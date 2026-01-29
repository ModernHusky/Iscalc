"""技能加载器模块

基于Agent Skills架构，从文件系统加载技能定义。
采用两阶段加载：
1. 启动阶段：只读取SKILL.md的YAML Frontmatter（元数据）
2. 执行阶段：按需读取完整SKILL.md内容
"""

import os
import re
import yaml
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any
from pathlib import Path


@dataclass
class SkillMetadata:
    """技能元数据（第一层 - 始终加载）"""
    name: str
    description: str
    keywords: List[str] = field(default_factory=list)
    applicable_types: List[str] = field(default_factory=list)
    match_rules: List[str] = field(default_factory=list)  # 正则表达式匹配规则
    path: str = ""  # SKILL.md文件路径


@dataclass
class SkillContent:
    """技能完整内容（第二层 - 按需加载）"""
    metadata: SkillMetadata
    full_content: str  # SKILL.md的Markdown正文


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
        self.skill_paths = self._determine_skill_paths(custom_skills_dir)
        self._skill_cache: Dict[str, SkillMetadata] = {}
        self._content_cache: Dict[str, SkillContent] = {}

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
    
    def discover_skills(self) -> List[SkillMetadata]:
        """发现所有技能（只读取元数据）
        
        扫描所有配置的路径，支持嵌套目录结构：
        - skills/rewrite/SKILL.md  (旧格式)
        - skills/commands/rewrite/SKILL.md  (新格式)
        
        如果同名技能出现在多个路径中，后扫描的（优先级高的）将覆盖先扫描的。
        """
        # 清空缓存以重新发现
        self._skill_cache = {}
        
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
        for skills_dir in self.skill_paths:
            scan_directory(skills_dir)
        
        return list(self._skill_cache.values())
    
    def _parse_frontmatter(self, skill_file: str) -> Optional[SkillMetadata]:
        """解析SKILL.md的YAML Frontmatter"""
        try:
            with open(skill_file, 'r', encoding='utf-8-sig') as f:
                content = f.read()
            # 规范化行尾 (CRLF -> LF)
            content = content.replace('\r\n', '\n').replace('\r', '\n')
            
            # 提取YAML Frontmatter
            match = re.match(r'^---\s*\n(.*?)\n---\s*\n', content, re.DOTALL)
            if not match:
                return None
            
            yaml_content = match.group(1)
            data = yaml.safe_load(yaml_content)
            
            return SkillMetadata(
                name=data.get('name', ''),
                description=data.get('description', ''),
                keywords=data.get('keywords', []),
                applicable_types=data.get('applicable_types', []),
                match_rules=data.get('match_rules', [])
            )
        except Exception as e:
            print(f"Warning: Failed to parse {skill_file}: {e}")
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
            return None
        
        try:
            with open(metadata.path, 'r', encoding='utf-8-sig') as f:
                content = f.read()
            # 规范化行尾 (CRLF -> LF)
            content = content.replace('\r\n', '\n').replace('\r', '\n')
            
            # 提取Markdown正文（去除Frontmatter）
            match = re.match(r'^---\s*\n.*?\n---\s*\n(.*)$', content, re.DOTALL)
            body = match.group(1) if match else content
            
            # 去除首尾空白
            body = body.strip()
            
            skill_content = SkillContent(metadata=metadata, full_content=body)
            self._content_cache[skill_name] = skill_content
            return skill_content
        except Exception as e:
            print(f"Warning: Failed to load {skill_name}: {e}")
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
            # 排除策略类技能，它们通常在第三层加载，不作为命令展示（或者可以展示，看设计）
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
                        flags = 0
                        # 如果配置是字典形式（新格式），提取 regex 和 flags
                        if isinstance(rule_config, dict):
                            pattern = rule_config.get('regex', '')
                            # 处理 flags (简单起见这里暂不处理 flags 字符串转换，通常 re.IGNORECASE 硬编码在下面)
                            # 如果需要支持 flags 列表，需解析字符串如 "re.DOTALL"
                        
                        if not isinstance(pattern, str):
                            continue

                        if re.search(pattern, text_to_match, re.IGNORECASE):
                            is_match = True
                            break
                    except re.error as e:
                        print(f"Warning: Invalid regex in skill {skill.name}: {pattern}, error: {e}")
            
            # 2. 如果没有 match_rules，则视为普通技能，不自动激活（或者可以添加其他逻辑）
            # 目前策略是：只有定义了 match_rules 的才会基于表达式自动激活
            # 对于 'strategy-examples' 这种，可能不需要基于表达式激活，而是由 prompts.py 显式调用
            
            if is_match:
                relevant.append(skill)
        
        return relevant
    
    def get_skill_instructions(self, skill_names: List[str]) -> str:
        """获取指定技能的完整指令
        
        按需加载第二层内容。
        """
        lines = []
        for name in skill_names:
            content = self.load_skill_content(name)
            if content:
                # 区分是普通技能还是策略指南
                # 策略指南通常不需要 "## 技能: name" 这样的标题，直接嵌入内容可能更自然
                # 但为了统一，我们保留简单的标题
                lines.append(f"\n## 技能: {name}\n")
                lines.append(content.full_content)
                lines.append("\n---\n")
        
        return "\n".join(lines)


# 全局技能加载器实例
_skill_loader: Optional[SkillLoader] = None


def get_skill_loader() -> SkillLoader:
    """获取全局技能加载器"""
    global _skill_loader
    if _skill_loader is None:
        _skill_loader = SkillLoader()
    return _skill_loader


def get_all_skill_metadata() -> str:
    """获取所有技能的元数据摘要"""
    return get_skill_loader().get_skills_summary()


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


# 在模块导入时自动加载
_ensure_skills_loaded()
