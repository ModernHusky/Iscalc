"""技能加载器模块

基于Agent Skills架构，从文件系统加载技能定义。
采用两阶段加载：
1. 启动阶段：只读取SKILL.md的YAML Frontmatter（元数据）
2. 执行阶段：按需读取完整SKILL.md内容

支持详细的日志调试输出。
"""

import os
import re
import yaml
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any, Tuple
from pathlib import Path

from .logger_config import get_phase_logger


@dataclass
class SkillMetadata:
    """技能元数据（第一层 - 始终加载）
    
    仅包含技能的基本信息，用于快速匹配和展示。
    Token 开销约 20-50。
    """
    name: str
    description: str
    keywords: List[str] = field(default_factory=list)
    applicable_types: List[str] = field(default_factory=list)
    match_rules: List[str] = field(default_factory=list)  # 正则表达式匹配规则
    path: str = ""  # SKILL.md 文件路径


@dataclass
class SkillContent:
    """技能完整内容（第二层 - 按需加载）
    
    包含 SKILL.md 的完整正文，提供详细指令和示例。
    Token 开销约 200-1000。
    """
    metadata: SkillMetadata
    full_content: str






class SkillLoader:
    """技能加载器
    
    实现渐进式披露：
    - discover_skills(): 扫描不同层级的技能目录（Project > Personal > Bundled）
    - load_skill_content(): 按需读取完整技能内容
    """
    
    def __init__(self, custom_skills_dir: Optional[str] = None):
        """初始化技能加载器
        
        Args:
            custom_skills_dir: 可选的自定义技能目录（用于测试或特定用途）
        """
        self.logger = get_phase_logger(__name__)
        self.skill_paths = self._determine_skill_paths(custom_skills_dir)

        # Skill caches
        self._skill_cache: Dict[str, SkillMetadata] = {}
        self._content_cache: Dict[str, SkillContent] = {}

        # Discovery cache (avoid expensive re-scan on every prompt build)
        self._discovery_done: bool = False
        self._last_discover_monotonic: float = 0.0
        # TTL is a pragmatic trade-off: avoids repeated filesystem walks while still allowing edits to take effect.
        self.discover_ttl_seconds: float = 30.0

        # Render caches (derived from metadata)
        self._skills_xml_cache: Optional[str] = None
        self._skills_categorized_xml_cache: Optional[str] = None

    def _determine_skill_paths(self, custom_dir: Optional[str]) -> List[str]:
        """确定技能扫描路径，按优先级从低到高排列（后加载覆盖先加载）
        
        优先级策略 (Project > Personal > Bundled):
        1. Bundled (Built-in): llm_iscalc/skills/
        2. Personal: ~/.claude/skills/
        3. Project: <project_root>/.claude/skills/
        """
        paths = []
        
        # 1. Bundled Skills (插件/内置层)
        # 默认为当前文件目录下的 skills 子目录
        bundled_path = os.path.join(os.path.dirname(__file__), "skills")
        if os.path.exists(bundled_path):
            paths.append(bundled_path)
        
        # 2. Personal Skills (用户层)
        try:
            home_dir = Path.home()
            personal_path = home_dir / ".claude" / "skills"
            if personal_path.exists():
                paths.append(str(personal_path))
        except Exception:
            pass  # 忽略获取 home 目录失败的情况

        # 3. Project Skills (项目层) - 优先级最高
        # 尝试寻找项目根目录 (假设含有 .git 或 .claude 目录，或者向上查找)
        project_root = self._find_project_root()
        if project_root:
            project_skill_path = os.path.join(project_root, ".claude", "skills")
            if os.path.exists(project_skill_path):
                paths.append(project_skill_path)
            
        # 如果提供了自定义目录，它具有最高优先级（通常用于测试）
        if custom_dir and os.path.exists(custom_dir):
            paths.append(custom_dir)
            
        return paths

    def _find_project_root(self) -> Optional[str]:
        """简单的项目根目录查找逻辑"""
        # 从当前文件所在目录开始向上查找
        current_dir = os.path.dirname(os.path.abspath(__file__))
        
        # 假设最多向上查5层
        for _ in range(5):
            # 检查特征文件/目录
            if os.path.exists(os.path.join(current_dir, ".git")) or \
               os.path.exists(os.path.join(current_dir, ".claude")):
                return current_dir
            
            parent_dir = os.path.dirname(current_dir)
            if parent_dir == current_dir: # 到达根目录
                break
            current_dir = parent_dir
            
        # Fallback: 如果这是作为包安装的，可能无法准确找到"项目根目录"
        # 此时可以假设运行时的 CWD 是项目根目录，或者返回 None
        # 这里为了安全，尝试返回 os.getcwd() 如果它看起来像个项目
        cwd = os.getcwd()
        if os.path.exists(os.path.join(cwd, ".git")):
            return cwd
            
        return None
    
    def discover_skills(self, force: bool = False) -> List[SkillMetadata]:
        """发现所有技能（只读取元数据）
        
        扫描所有配置的路径，支持嵌套目录结构：
        - skills/rewrite/SKILL.md  (旧格式)
        - skills/commands/rewrite/SKILL.md  (新格式)
        
        如果同名技能出现在多个路径中，后扫描的（优先级高的）将覆盖先扫描的。
        
        Notes:
            此操作会进行文件系统遍历与 YAML frontmatter 解析，比较昂贵。
            默认启用 TTL 缓存，避免在每次 LLM prompt build 时都重复扫描。
        """
        import time

        now = time.monotonic()
        if not force and self._discovery_done and (now - self._last_discover_monotonic) < self.discover_ttl_seconds:
            return list(self._skill_cache.values())

        # Rescan: clear caches
        self._skill_cache = {}
        self._content_cache = {}
        self._skills_xml_cache = None
        self._skills_categorized_xml_cache = None
        
        def scan_directory(base_dir: str, depth: int = 0):
            """递归扫描目录，最多两层嵌套"""
            if not os.path.exists(base_dir) or depth > 2:
                return
                
            for item_name in os.listdir(base_dir):
                item_path = os.path.join(base_dir, item_name)
                
                if not os.path.isdir(item_path):
                    continue
                    
                skill_file = os.path.join(item_path, "SKILL.md")
                
                if os.path.exists(skill_file):
                    # 找到 SKILL.md，解析元数据
                    metadata = self._parse_frontmatter(skill_file)
                    if metadata:
                        metadata.path = skill_file
                        self._skill_cache[metadata.name] = metadata
                else:
                    # 没有 SKILL.md，可能是分类目录，继续递归
                    scan_directory(item_path, depth + 1)
        
        # 按照 paths 列表顺序扫描（Bundled -> Personal -> Project）
        self.logger.discovery("开始扫描技能目录...")
        
        path_labels = ["Bundled", "Personal", "Project", "Custom"]
        for i, skills_dir in enumerate(self.skill_paths):
            label = path_labels[i] if i < len(path_labels) else f"Path-{i}"
            self.logger.info("   └── [%s] %s", label, skills_dir)
            scan_directory(skills_dir)
        
        skills = list(self._skill_cache.values())
        self.logger.discovery("发现 %d 个技能", len(skills))
        
        # 按类别分组输出
        if skills:
            categories = {}
            for s in skills:
                # 从路径提取类别
                parts = s.path.replace("\\", "/").split("/")
                if "commands" in parts:
                    cat = "commands"
                elif "states" in parts:
                    cat = "states"
                elif "strategies" in parts:
                    cat = "strategies"
                elif "general" in parts:
                    cat = "general"
                else:
                    cat = "other"
                categories.setdefault(cat, []).append(s.name)
            
            for cat, names in categories.items():
                self.logger.info("   └── [%s]: %s", cat, ", ".join(names[:5]) + ("..." if len(names) > 5 else ""))

        # Mark discovery as fresh
        self._discovery_done = True
        self._last_discover_monotonic = now

        return skills
    
    def _parse_frontmatter(self, skill_file: str) -> Optional[SkillMetadata]:
        """解析SKILL.md的YAML Frontmatter (高效读取)"""
        try:
            yaml_lines = []
            with open(skill_file, 'r', encoding='utf-8-sig') as f:
                # 检查第一行是否是 ---
                first_line = f.readline()
                if not first_line.strip() == '---':
                    return None
                
                # 读取直到下一个 ---
                for line in f:
                    if line.strip() == '---':
                        break
                    yaml_lines.append(line)
            
            if not yaml_lines:
                return None
                
            yaml_content = ''.join(yaml_lines)
            data = yaml.safe_load(yaml_content)
            
            return SkillMetadata(
                name=data.get('name', ''),
                description=data.get('description', ''),
                keywords=data.get('keywords', []),
                applicable_types=data.get('applicable_types', []),
                match_rules=data.get('match_rules', [])
            )
        except Exception as e:
            self.logger.warning("解析技能文件失败: %s - %s", skill_file, str(e))
            return None
    
    def load_skill_content(self, skill_name: str) -> Optional[SkillContent]:
        """加载技能完整内容（按需）
        
        当需要使用某个技能时，读取完整的SKILL.md内容。
        """
        if skill_name in self._content_cache:
            return self._content_cache[skill_name]
        
        if skill_name not in self._skill_cache:
            self.discover_skills()
        
        metadata = self._skill_cache.get(skill_name)
        if not metadata or not metadata.path:
            self.logger.warning(f"⚠️ Load Content: Skill '{skill_name}' not found or no path.")
            return None
        
        try:
            # 技能内容加载日志
            self.logger.loading("加载技能内容: %s", skill_name)
            self.logger.info("   路径: %s", metadata.path)
            
            content_lines = []
            with open(metadata.path, 'r', encoding='utf-8-sig') as f:
                # 跳过 Frontmatter
                first_line = f.readline()
                if first_line.strip() == '---':
                    for line in f:
                        if line.strip() == '---':
                            break
                
                # 读取剩余内容
                content_lines = f.readlines()
                
            body = "".join(content_lines).strip()
            
            skill_content = SkillContent(metadata=metadata, full_content=body)
            self._content_cache[skill_name] = skill_content
            self.logger.skill_load("技能内容已加载: %s (%d 字节)", skill_name, len(body))
            return skill_content
        except Exception as e:
            self.logger.error_phase("加载技能内容失败: %s - %s", skill_name, str(e))
            return None
    

    def get_skills_summary(self) -> str:
        """获取所有技能的摘要（第一层）
        
        用于注入系统提示。
        """
        skills = self.discover_skills()
        if not skills:
            return ""
        
        lines = ["## 可用技能\n"]
        for skill in skills:
            # 这里我们选择只展示描述，让LLM知道有什么能力
            lines.append(f"- **{skill.name}**: {skill.description}")
        
        return "\n".join(lines)
    
    def get_relevant_skills(self, expression: str, user_instruction: Optional[str] = None) -> List[SkillMetadata]:
        """根据表达式和用户指令获取相关技能
        
        使用match_rules正则匹配。
        """
        skills = self.discover_skills()
        relevant = []
        
        # 合并匹配文本，加入换行符分隔以避免意外拼接
        text_to_match = expression
        if user_instruction:
            text_to_match += "\n" + user_instruction

        for skill in skills:
            is_match = False
            # 1. 检查 match_rules
            if skill.match_rules:
                for rule_config in skill.match_rules:
                    try:
                        pattern = rule_config
                        if isinstance(rule_config, dict):
                            pattern = rule_config.get('regex', '')
                        
                        if not isinstance(pattern, str):
                            continue

                        if re.search(pattern, text_to_match, re.IGNORECASE):
                            is_match = True
                            break
                    except re.error as e:
                        print(f"Warning: Invalid regex in skill {skill.name}: {rule_config}, error: {e}")
            
            # 2. 如果没有 match_rules，则视为普通技能，不自动激活（或者可以添加其他逻辑）
            # 目前策略是：只有定义了 match_rules 的才会基于表达式自动激活
            # 对于 'strategy-examples' 这种，可能不需要基于表达式激活，而是由 prompts.py 显式调用
            
            if is_match:
                relevant.append(skill)
        
        return relevant
    
    def search_skills(self, query: str, limit: int = 5) -> List[SkillMetadata]:
        """在技能元数据上做轻量关键词/模糊搜索（不加载全文）。

        用途：Search-o1 风格的 `<|load_skill|>...<|end_load_skill|>` 标记允许写“关键词”，
        系统可据此选择最相关的技能文件。
        """
        import re

        q = (query or "").strip()
        if not q:
            return []

        # Strip common prefixes (model may output `search: xxx`)
        q_lower = q.lower()
        for prefix in ("search:", "query:", "kw:", "skill:", "关键词:", "关键字:"):
            if q_lower.startswith(prefix):
                q = q[len(prefix):].strip()
                q_lower = q.lower()
                break

        skills = self.discover_skills()

        # Tokenize conservatively (keep hyphenated words meaningful)
        tokens = [t for t in re.split(r"[\s_/]+|:+", q_lower) if t]
        if not tokens:
            tokens = [q_lower]

        scored: List[Tuple[float, SkillMetadata]] = []
        for s in skills:
            name = (s.name or "").lower()
            desc = (s.description or "").lower()
            kws = [k.lower() for k in (s.keywords or []) if isinstance(k, str)]

            score = 0.0

            # Strong signals
            if name == q_lower:
                score += 100.0
            if q_lower.replace(" ", "-") == name:
                score += 80.0
            if q_lower in name:
                score += 30.0

            # Token matches
            for t in tokens:
                if t in name:
                    score += 10.0
                if t in desc:
                    score += 2.0
                if any(t in kw for kw in kws):
                    score += 6.0

            # Bonus: regex match_rules can hint intent (best-effort)
            for rule_config in (s.match_rules or []):
                try:
                    pattern = rule_config.get('regex', '') if isinstance(rule_config, dict) else rule_config
                    if isinstance(pattern, str) and pattern and re.search(pattern, query, re.IGNORECASE):
                        score += 1.0
                        break
                except Exception:
                    continue

            if score > 0.0:
                scored.append((score, s))

        scored.sort(key=lambda x: x[0], reverse=True)
        return [s for _, s in scored[:limit]]
    
    
    def get_skill_instructions(self, skill_names: List[str]) -> str:
        """获取指定技能的完整指令
        
        按需加载第二层内容。
        """
        lines = []
        for name in skill_names:
            content = self.load_skill_content(name)
            if content:
                lines.append(f"\n## 技能: {name}\n")
                lines.append(content.full_content)
                lines.append("\n---\n")
        
        return "\n".join(lines)
    
    def get_skills_xml(self) -> str:
        """获取所有技能的XML格式元数据（第一层）"""
        skills = self.discover_skills()
        if self._skills_xml_cache is not None:
            return self._skills_xml_cache

        if not skills:
            self._skills_xml_cache = "<skill_list></skill_list>"
            return self._skills_xml_cache
        
        # Calculate base directory for relative paths
        base_dir = os.path.dirname(os.path.abspath(__file__))
        
        xml_lines = ["<skill_list>"]
        for skill in skills:
            # XML Escape for description
            desc = skill.description.replace("<", "&lt;").replace(">", "&gt;")
            # Calculate relative path
            try:
                rel_path = os.path.relpath(skill.path, base_dir).replace("\\", "/")
            except ValueError:
                rel_path = skill.path # Fallback to absolute if on different drive
                
            xml_lines.append(f'    <skill name="{skill.name}" path="{rel_path}">{desc}</skill>')
        xml_lines.append("</skill_list>")
        self._skills_xml_cache = "\n".join(xml_lines)
        return self._skills_xml_cache
    
    def get_skill_path(self, skill_name: str) -> Optional[str]:
        """获取技能文件的绝对路径"""
        if skill_name not in self._skill_cache:
            self.discover_skills()
        
        meta = self._skill_cache.get(skill_name)
        return meta.path if meta else None


# 全局技能加载器实例
_skill_loader: Optional[SkillLoader] = None


def get_skill_loader() -> SkillLoader:
    """获取全局技能加载器"""
    global _skill_loader
    if _skill_loader is None:
        _skill_loader = SkillLoader()
    return _skill_loader


def predict_skills_for_expression(expression: str, state: str = "CALCULATE") -> List[str]:
    """根据表达式特征和状态预测需要的技能
    
    Args:
        expression: 数学表达式
        state: 当前状态 (CALCULATE, PROVE, INDUCTION 等)
    
    Returns:
        建议预加载的技能名称列表
    """
    skills = []
    expr_lower = expression.lower()
    
    # 状态相关技能
    state_skill_map = {
        "PROVE": ["state-prove"],
        "INDUCTION": ["state-induction"],
        "CALCULATE": ["state-calculate"]
    }
    if state in state_skill_map:
        skills.extend(state_skill_map[state])
    
    # 表达式模式匹配
    if "int " in expr_lower or "∫" in expression:
        skills.append("strategy-integral")
        if any(x in expr_lower for x in ["sin", "cos", "tan"]):
            skills.append("integrate-by-parts")
        if "/" in expression or "1/(" in expr_lower:
            skills.append("partial-fraction")
    
    if "lim" in expr_lower or "→" in expression:
        skills.append("strategy-limit")
    
    if any(x in expr_lower for x in ["sin", "cos", "tan", "arcsin", "arccos"]):
        skills.append("rewrite")
    
    if "sum(" in expr_lower or "∑" in expression:
        skills.append("strategy-summation")
    
    # 去重并限制数量
    seen = set()
    unique_skills = []
    for s in skills:
        if s not in seen:
            seen.add(s)
            unique_skills.append(s)
    
    return unique_skills[:3]  # 最多预加载3个技能


def get_all_skill_metadata() -> str:
    """获取所有技能的元数据摘要（已废弃，建议使用 get_skills_xml）"""
    return get_skill_loader().get_skills_summary()


def get_skills_xml() -> str:
    """获取所有技能的XML元数据"""
    return get_skill_loader().get_skills_xml()


def get_skills_categorized_xml() -> str:
    """获取按类别分组的技能列表（增强版）
    
    返回更易于 LLM 理解的分类格式。
    """
    loader = get_skill_loader()
    skills = loader.discover_skills()

    if getattr(loader, "_skills_categorized_xml_cache", None) is not None:
        return loader._skills_categorized_xml_cache  # type: ignore[attr-defined]

    if not skills:
        loader._skills_categorized_xml_cache = "<skill_categories></skill_categories>"  # type: ignore[attr-defined]
        return loader._skills_categorized_xml_cache  # type: ignore[attr-defined]
    
    # Calculate base directory for relative paths
    base_dir = os.path.dirname(os.path.abspath(__file__))
    
    # 按类别分组
    categories = {
        "strategies": {"name": "策略技能", "desc": "处理特定类型问题的整体策略，建议首先阅读", "skills": []},
        "commands": {"name": "命令技能", "desc": "具体的 iscalc 命令用法", "skills": []},
        "states": {"name": "状态技能", "desc": "不同求解状态下的操作指南", "skills": []},
        "general": {"name": "通用技能", "desc": "通用工具和示例", "skills": []},
    }
    
    for skill in skills:
        # 从路径提取类别
        parts = skill.path.replace("\\", "/").split("/")
        if "strategies" in parts:
            cat = "strategies"
        elif "commands" in parts:
            cat = "commands"
        elif "states" in parts:
            cat = "states"
        elif "general" in parts:
            cat = "general"
        else:
            cat = "general"
        
        # Calculate relative path
        try:
            rel_path = os.path.relpath(skill.path, base_dir).replace("\\", "/")
        except ValueError:
            rel_path = skill.path
        
        # XML Escape
        desc = skill.description.replace("<", "&lt;").replace(">", "&gt;")
        
        categories[cat]["skills"].append({
            "name": skill.name,
            "path": rel_path,
            "desc": desc
        })
    
    # 生成 XML
    xml_lines = ["<skill_categories>"]
    
    # 按优先级顺序输出：策略 > 命令 > 状态 > 通用
    for cat_key in ["strategies", "commands", "states", "general"]:
        cat = categories[cat_key]
        if cat["skills"]:
            xml_lines.append(f'  <category name="{cat["name"]}" description="{cat["desc"]}">')
            for s in cat["skills"]:
                xml_lines.append(f'    <skill name="{s["name"]}" path="{s["path"]}">{s["desc"]}</skill>')
            xml_lines.append("  </category>")
    
    xml_lines.append("</skill_categories>")
    loader._skills_categorized_xml_cache = "\n".join(xml_lines)  # type: ignore[attr-defined]
    return loader._skills_categorized_xml_cache  # type: ignore[attr-defined]


def get_relevant_skills(expression: str, user_instruction: Optional[str] = None) -> List[SkillMetadata]:
    """根据表达式和用户指令获取相关技能"""
    return get_skill_loader().get_relevant_skills(expression, user_instruction)


def get_skill_details(skills: List[SkillMetadata], include_examples: bool = False) -> str:
    """获取技能详情"""
    loader = get_skill_loader()
    skill_names = [s.name for s in skills]
    return loader.get_skill_instructions(skill_names)



def get_state_skills(state_name: str) -> List[SkillMetadata]:
    """根据状态名称获取相关的状态技能
    
    Args:
        state_name: 状态名称，如 'CALCULATE', 'PROVE', 'INDUCTION' 等（大小写不敏感）
    
    Returns:
        匹配的状态技能列表
    """
    loader = get_skill_loader()
    all_skills = loader.discover_skills()
    
    state_name_lower = state_name.lower()
    relevant = []
    
    # 状态名称到技能名称的映射
    state_skill_map = {
        'calculate': ['state-calculate', 'state-done'],
        'prove': ['state-prove', 'state-done', 'proof-states'],  # 添加state-prove
        'induction': ['state-induction', 'state-done', 'induction-proof'],
        'case': ['state-done', 'case-analysis'],
        'initial': ['proof-states', 'state-calculate'],
    }
    
    # 获取该状态应该加载的技能名称
    target_skill_names = state_skill_map.get(state_name_lower, [])
    
    for skill in all_skills:
        # 1. 检查技能名是否在映射列表中
        if skill.name in target_skill_names:
            relevant.append(skill)
        # 2. 检查技能名是否包含状态名（如 state-induction）
        elif state_name_lower in skill.name.lower():
            if skill not in relevant:
                relevant.append(skill)
    
    return relevant


# 为了向后兼容，保留 COMMAND_SKILLS
# 这些将从文件系统动态加载
@dataclass
class CommandSkill:
    """命令技能定义（兼容旧接口）"""
    name: str
    summary: str
    keywords: List[str] = field(default_factory=list)
    applicable_types: List[str] = field(default_factory=list)
    
    def get_metadata(self) -> str:
        return f"- **{self.name}**: {self.summary}"
    
    def get_core_instruction(self) -> str:
        return f"### {self.name}\n{self.summary}"
    
    def get_full_detail(self) -> str:
        return self.get_core_instruction()


def _load_command_skills() -> List[CommandSkill]:
    """从文件系统加载技能，转换为CommandSkill格式"""
    loader = get_skill_loader()
    metadata_list = loader.discover_skills()
    
    skills = []
    for meta in metadata_list:
        skill = CommandSkill(
            name=meta.name,
            summary=meta.description,
            keywords=meta.keywords,
            applicable_types=meta.applicable_types
        )
        skills.append(skill)
    
    return skills


# 惰性加载
COMMAND_SKILLS: List[CommandSkill] = []


def _ensure_skills_loaded():
    global COMMAND_SKILLS
    if not COMMAND_SKILLS:
        COMMAND_SKILLS = _load_command_skills()


# 注意：不要在模块导入时自动扫描技能目录。
# 如需旧接口，请显式调用 _ensure_skills_loaded()。


def detect_skill_mentions(text: str, loaded_skills: Optional[List[str]] = None) -> List[Dict[str, str]]:
    """从文本中检测提到的技能名称
    
    扫描 LLM 的 thinking 内容，检测是否提到了某个技能名称。
    如果提到了，返回该技能的相关信息以便自动加载。
    
    Args:
        text: LLM 的 thinking 或其他文本内容
        loaded_skills: 已加载的技能路径列表，用于避免重复加载
        
    Returns:
        包含检测到的技能信息的列表 [{"name": "xxx", "path": "...", "full_path": "..."}]
    """
    if not text:
        return []
    
    loaded_skills = loaded_skills or []
    detected = []
    
    loader = get_skill_loader()
    all_skills = loader.discover_skills()
    
    # 构建技能名称到元数据的映射
    skill_map = {}
    for skill in all_skills:
        # 使用多种可能的匹配模式
        # 1. 技能名称本身 (如 "merge-evalat")
        skill_map[skill.name.lower()] = skill
        # 2. 文件路径的变体 (如 "skills/merge-evalat")
        rel_path = os.path.relpath(skill.path, os.path.dirname(os.path.dirname(skill.path)))
        skill_map[rel_path.replace("\\", "/").lower()] = skill
        # 3. 简短变体 (如 "merge evalat", "mergeevalat")
        skill_map[skill.name.replace("-", " ").lower()] = skill
        skill_map[skill.name.replace("-", "").lower()] = skill
    
    text_lower = text.lower()
    
    # 检测匹配
    for pattern, skill in skill_map.items():
        if pattern in text_lower:
            # 检查是否已加载
            if skill.path in loaded_skills:
                continue
            
            # 检查是否已在检测结果中
            if any(d["path"] == skill.path for d in detected):
                continue
            
            # 计算相对路径（相对于 llm_iscalc 目录）
            base_dir = os.path.dirname(os.path.abspath(__file__))
            try:
                rel_path = os.path.relpath(skill.path, base_dir).replace("\\", "/")
            except ValueError:
                rel_path = skill.path
            
            detected.append({
                "name": skill.name,
                "path": rel_path,  # 相对路径，用于 cat 命令
                "full_path": skill.path  # 绝对路径，用于实际加载
            })
    
    return detected
