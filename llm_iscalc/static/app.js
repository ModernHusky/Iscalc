// 调试信息
console.log('app.js 已加载');

// 全局状态
let isSolving = false;
let eventSource = null;
let loadedSkillsThisSession = new Set(); // 记录本次会话已通知的技能

// ============ Search-o1 风格：技能加载通知 ============

/**
 * 显示技能加载通知（右上角滑入）
 * @param {string} skillName - 技能名称
 */
function showSkillLoadNotification(skillName) {
    const container = document.getElementById('skill-notifications');
    if (!container) return;

    // 避免重复通知
    if (loadedSkillsThisSession.has(skillName)) return;
    loadedSkillsThisSession.add(skillName);

    // 创建通知元素
    const notification = document.createElement('div');
    notification.className = 'skill-load-notification';
    notification.innerHTML = `
        <span class="icon">🔖</span>
        <span>已加载技能: <strong>${skillName}</strong></span>
    `;

    container.appendChild(notification);

    // 3秒后淡出
    setTimeout(() => {
        notification.classList.add('fade-out');
        setTimeout(() => {
            notification.remove();
        }, 300);
    }, 3000);
}

/**
 * 从文本中检测技能加载标记
 * @param {string} text - LLM输出的文本
 */
function detectSkillLoading(text) {
    const loadPattern = /\[✓ 已加载技能: (.+?)\]/g;
    const matches = [...text.matchAll(loadPattern)];
    matches.forEach(match => {
        const skillName = match[1];
        showSkillLoadNotification(skillName);
    });
}

// DOM 元素
const expressionInput = document.getElementById('expression-input');
const conditionsInput = document.getElementById('conditions-input');
const instructionInput = document.getElementById('instruction-input');
const solveBtn = document.getElementById('solve-btn');
const stopBtn = document.getElementById('stop-btn');
const clearBtn = document.getElementById('clear-btn');
const solveBtnText = document.getElementById('solve-btn-text');
const solveSpinner = document.getElementById('solve-spinner');

console.log('DOM 元素获取:', {
    expressionInput: !!expressionInput,
    conditionsInput: !!conditionsInput,
    solveBtn: !!solveBtn,
    stopBtn: !!stopBtn,
    clearBtn: !!clearBtn
});

const resultOutput = document.getElementById('result-output');
const thinkingOutput = document.getElementById('thinking-output');
const commandOutput = document.getElementById('command-output');
const errorOutput = document.getElementById('error-output');
const stdCommandList = document.getElementById('std-command-list');
const clearThinkingBtn = document.getElementById('clear-thinking-btn');
const historyList = document.getElementById('history-list');
const clearHistoryBtn = document.getElementById('clear-history-btn');

// 历史记录管理
let solveHistory = [];

// 从localStorage加载历史记录
function loadHistory() {
    try {
        const saved = localStorage.getItem('llm_iscalc_history');
        if (saved) {
            solveHistory = JSON.parse(saved);
            renderHistory();
        }
    } catch (e) {
        console.error('加载历史记录失败:', e);
    }
}

// 保存历史记录到localStorage
function saveHistory() {
    try {
        // 只保留最近20条
        if (solveHistory.length > 20) {
            solveHistory = solveHistory.slice(-20);
        }
        localStorage.setItem('llm_iscalc_history', JSON.stringify(solveHistory));
    } catch (e) {
        console.error('保存历史记录失败:', e);
    }
}

// 创建历史记录项 DOM
function createHistoryItemDOM(item) {
    const div = document.createElement('div');
    div.className = 'history-item';
    // 设置唯一 ID (用于动画追踪)
    div.dataset.key = `${item.expression}||${item.conditions || ''}`;

    div.innerHTML = `
        <div class="history-item-expr" title="${item.expression}">${item.expression}</div>
        ${item.conditions ? `<div class="history-item-cond" title="${item.conditions}">条件: ${item.conditions}</div>` : ''}
        <div class="history-item-time">${item.timestamp}</div>
    `;

    div.addEventListener('click', () => {
        expressionInput.value = item.expression;
        conditionsInput.value = item.conditions;
        expressionInput.focus();
    });

    return div;
}

// 添加历史记录 (带动画)
function addHistory(expression, conditions) {
    conditions = conditions || '';

    const item = {
        expression: expression,
        conditions: conditions,
        timestamp: new Date().toLocaleString('zh-CN', {
            month: '2-digit',
            day: '2-digit',
            hour: '2-digit',
            minute: '2-digit'
        })
    };

    // 1. 记录初始位置 (FLIP - First)
    const firstPositions = new Map();
    const items = historyList.querySelectorAll('.history-item');
    items.forEach(el => {
        if (el.dataset.key) {
            firstPositions.set(el.dataset.key, el.getBoundingClientRect());
        }
    });

    // 2. 更新数据模型
    const key = `${expression}||${conditions}`;
    const existingIndex = solveHistory.findIndex(h =>
        h.expression === expression && h.conditions === conditions
    );

    if (existingIndex !== -1) {
        solveHistory.splice(existingIndex, 1);
    }
    solveHistory.push(item);

    // 确保不超过限制
    if (solveHistory.length > 20) {
        solveHistory.shift(); // 移除最早的 (注意: 渲染是倒序, 所以数据头部是最早的)
    }
    saveHistory();

    // 3. 更新 DOM (FLIP - Last)
    // 移除空状态
    const emptyState = historyList.querySelector('.text-gray-400');
    if (emptyState) emptyState.remove();

    // 查找已存在的 DOM 元素
    let domElement = null;
    items.forEach(el => {
        if (el.dataset.key === key) domElement = el;
    });

    if (domElement) {
        // 移动现有元素到顶部 (prepend)
        // 更新时间
        const timeEl = domElement.querySelector('.history-item-time');
        if (timeEl) timeEl.textContent = item.timestamp;

        historyList.prepend(domElement);
    } else {
        // 创建新元素
        domElement = createHistoryItemDOM(item);
        historyList.prepend(domElement);
    }

    // 移除多余的 DOM 元素
    while (historyList.children.length > 20) {
        historyList.lastElementChild.remove();
    }

    // 4. 执行动画 (FLIP - Invert & Play)
    requestAnimationFrame(() => {
        const allItems = historyList.querySelectorAll('.history-item');

        allItems.forEach(el => {
            const firstRect = firstPositions.get(el.dataset.key);
            const lastRect = el.getBoundingClientRect();

            if (firstRect) {
                // 移动的元素
                const deltaY = firstRect.top - lastRect.top;

                // 只有当位置确实发生改变时才动画
                if (Math.abs(deltaY) > 0) {
                    el.style.transition = 'none';
                    el.style.transform = `translateY(${deltaY}px)`;

                    // Trigger reflow
                    el.offsetHeight;

                    // Play
                    el.style.transition = 'transform 0.4s cubic-bezier(0.25, 1, 0.5, 1)';
                    el.style.transform = '';
                }
            } else {
                // 新进入的元素 (Fade in + Slide down)
                el.style.transition = 'none';
                el.style.opacity = '0';
                el.style.transform = 'translateY(-20px)';

                // Trigger reflow
                el.offsetHeight;

                el.style.transition = 'all 0.4s cubic-bezier(0.25, 1, 0.5, 1)';
                el.style.opacity = '1';
                el.style.transform = '';
            }

            // 动画结束后清理 style
            const onEnd = () => {
                el.style.transition = '';
                el.style.transform = '';
                el.style.opacity = '';
                el.removeEventListener('transitionend', onEnd);
            };
            el.addEventListener('transitionend', onEnd);
        });
    });
}

// 渲染历史记录
function renderHistory() {
    if (solveHistory.length === 0) {
        historyList.innerHTML = '<div class="text-gray-400 text-center text-xs py-2">暂无历史记录</div>';
        return;
    }

    historyList.innerHTML = '';

    // 倒序显示（最新的在上面）
    // 注意：addHistory 已经改为增量更新 DOM，renderHistory 仅用于初始化加载
    for (let i = solveHistory.length - 1; i >= 0; i--) {
        const item = solveHistory[i];
        const div = createHistoryItemDOM(item);
        historyList.appendChild(div);
    }
}

// 清空历史记录
clearHistoryBtn.addEventListener('click', () => {
    if (confirm('确定要清空所有历史记录吗？')) {
        solveHistory = [];
        saveHistory();
        renderHistory();
    }
});

// 页面加载时加载历史记录
loadHistory();

// 智能滚动管理器
class SmartScrollManager {
    constructor(element) {
        this.element = element;
        this.userIsScrolling = false;
        this.scrollTimeout = null;
        this.attachListeners();
    }

    attachListeners() {
        this.element.addEventListener('scroll', () => {
            // 用户手动滚动时，标记为正在滚动
            this.userIsScrolling = true;

            // 清除之前的超时
            if (this.scrollTimeout) {
                clearTimeout(this.scrollTimeout);
            }

            // 如果用户滚动到底部，恢复自动滚动
            if (this.isScrolledToBottom()) {
                this.userIsScrolling = false;
            } else {
                // 2秒后如果没有新的滚动事件，恢复自动滚动
                this.scrollTimeout = setTimeout(() => {
                    if (this.isScrolledToBottom()) {
                        this.userIsScrolling = false;
                    }
                }, 2000);
            }
        });
    }

    isScrolledToBottom() {
        const threshold = 5; // 5px 的容差，减少抖动（原为50px）
        return this.element.scrollHeight - this.element.scrollTop - this.element.clientHeight < threshold;
    }

    scrollToBottom() {
        // 只有当用户没有正在滚动时才自动滚动
        if (!this.userIsScrolling) {
            this.element.scrollTop = this.element.scrollHeight;
        }
    }

    reset() {
        this.userIsScrolling = false;
        if (this.scrollTimeout) {
            clearTimeout(this.scrollTimeout);
            this.scrollTimeout = null;
        }
    }
}

// 初始化各个区域的滚动管理器
const scrollManagers = {
    thinking: new SmartScrollManager(thinkingOutput),
    result: new SmartScrollManager(resultOutput),
    command: new SmartScrollManager(commandOutput),
    error: new SmartScrollManager(errorOutput),
    stdCommand: new SmartScrollManager(stdCommandList) // 侧边栏也加上
};

// 统一的智能滚动调用接口
function smartScroll(element) {
    // 查找对应的管理器
    let manager = null;
    if (element === thinkingOutput) manager = scrollManagers.thinking;
    else if (element === resultOutput) manager = scrollManagers.result;
    else if (element === commandOutput) manager = scrollManagers.command;
    else if (element === errorOutput) manager = scrollManagers.error;
    else if (element === stdCommandList) manager = scrollManagers.stdCommand;

    if (manager) {
        manager.scrollToBottom();
    } else {
        // Fallback
        element.scrollTop = element.scrollHeight;
    }
}

let thinkingStepCount = 0;
let currentCard = null; // 当前正在更新的卡片
let currentSidebarItem = null; // Track current command item in sidebar
let currentFields = {}; // 当前卡片的各个字段元素
let lastThinkingValue = ''; // 记录上一次的 thinking 值
let hasContent = false; // 标记是否已有内容
let lastStep = -1; // 记录上一次的步骤号

// 清空思考历史
clearThinkingBtn.addEventListener('click', () => {
    thinkingStepCount = 0;
    currentCard = null;
    currentFields = {};
    lastThinkingValue = '';
    hasContent = false;
    lastStep = -1; // 重置步骤号
    scrollManagers.thinking.reset(); // 重置滚动状态
    thinkingOutput.innerHTML = '<div class="text-gray-400 text-center py-8 text-xs">AI 思考过程将在这里显示...</div>';
});

// 字段定义配置
const FIELD_CONFIG = {
    'thinking': { icon: '🧠', label: '分析' },
    'command': { icon: '⚡', label: '命令' },
    'explanation': { icon: '📝', label: '解释' },
    'is_final': { icon: '✓', label: '证明结束否' }
};

// 创建新的思考卡片
function createNewCard() {
    thinkingStepCount++;

    // 清除占位文本
    if (thinkingOutput.querySelector('.text-gray-400')) {
        thinkingOutput.innerHTML = '';
    }

    // 创建卡片容器
    const card = document.createElement('div');
    card.className = 'thinking-card';
    card.innerHTML = '<div class="thinking-card-grid"></div>';

    thinkingOutput.appendChild(card);

    // 重置当前卡片引用
    currentCard = card;
    currentSidebarItem = null; // Reset for new card
    currentFields = {}; // 清空字段引用

    // 注意：不再预先创建所有字段。字段将在 updateField 中按需创建。
    return card;
}

// 添加或更新字段
function updateField(fieldName, content) {
    if (!currentCard) return;

    // 如果字段尚未创建，则根据配置动态创建
    if (!currentFields[fieldName]) {
        const config = FIELD_CONFIG[fieldName];
        if (!config) return; // 未知字段，忽略

        const grid = currentCard.querySelector('.thinking-card-grid');
        const fieldDiv = document.createElement('div');
        fieldDiv.className = `thinking-field ${fieldName}`;
        // 初始状态可以添加简单的淡入动画效果
        fieldDiv.style.animation = 'fadeIn 0.3s ease-out';

        fieldDiv.innerHTML = `
            <div class="thinking-field-label">
                <span>${config.icon}</span>
                <span>${config.label}</span>
            </div>
            <div class="thinking-field-content"></div>
        `;
        grid.appendChild(fieldDiv);
        currentFields[fieldName] = fieldDiv.querySelector('.thinking-field-content');

        // 当新字段出现时，确保滚动到底部
        smartScroll(thinkingOutput);
    }

    // 更新内容
    const contentDiv = currentFields[fieldName];
    if (fieldName === 'is_final') {
        const isTrue = content === 'true';
        contentDiv.textContent = isTrue ? '✓ 是' : '✗ 否';
        if (isTrue) {
            contentDiv.parentElement.classList.add('is-true');
        }
    } else {
        // 流式更新：直接显示最新内容
        contentDiv.textContent = content.replace(/\\n/g, '\n');
    }

    // 智能滚动
    smartScroll(thinkingOutput);
}

// 从结构化对象更新思考内容（系统层面解决方案）
function updateThinkingFromObject(thinkingObj) {
    console.log('=== 更新思考对象 ===');
    console.log('完整对象:', thinkingObj);

    const thinking = thinkingObj.thinking || '';
    const command = thinkingObj.command || '';
    const explanation = thinkingObj.explanation || '';
    const isFinal = thinkingObj.is_final; // Don't default to false yet, to distinguish undefined
    const step = thinkingObj.step || 0;

    // 检查步骤号
    if (step < 1) return;

    const hasAnyContent = thinking || command || explanation || (isFinal !== undefined);

    console.log('提取字段:', { thinking: thinking.substring(0, 50), command, explanation: explanation.substring(0, 50), isFinal, step });
    console.log('当前状态:', { lastStep, currentCard: !!currentCard });

    // 判断是否需要创建新卡片
    // 如果步骤变化，尝试创建新卡片
    if (step !== lastStep) {
        // 只有当有实质内容时才创建卡片，避免空白卡片
        if (hasAnyContent) {
            createNewCard(); // createNewCard doesn't take step as argument
            lastStep = step;
            console.log('✓ 新卡片已创建，lastStep更新为:', lastStep);
        } else {
            console.log('→ 步骤变化但无实质内容，等待内容');
            return; // 等待内容
        }
    } else {
        console.log('→ 同一步骤，更新当前卡片');
    }

    // 在当前卡片内更新字段（流式更新）
    // 始终调用 updateField，即使内容为空，以确保字段存在（如果卡片已创建）
    // 但此时 updateField 会更新已存在的卡片
    if (currentCard) {
        if (thinking) {
            updateField('thinking', thinking);
            lastThinkingValue = thinking;
            hasContent = true;
        }
        if (command) {
            updateField('command', command);

            // Sync to sidebar
            // 所有命令都添加到sidebar,包括中间步骤
            // Filter out "None" or empty commands
            if (command && command !== "None" && command.trim() !== "") {
                // 注意:不再区分用户命令和中间步骤,所有命令都添加
                // 因为中间步骤现在会通过独立的COMMAND事件发送
                if (!currentSidebarItem) {
                    currentSidebarItem = addCommandToSidebar(command);
                } else {
                    updateCommandSidebarItem(currentSidebarItem, command);
                }
            }

            hasContent = true;
        }
        if (explanation) {
            updateField('explanation', explanation);
            hasContent = true;
        }
        // isFinal 是布尔值，需要特殊处理
        if (isFinal !== undefined) {
            updateField('is_final', String(isFinal));
            hasContent = true;

            // If final and no command yet for this step, add 'simplify' to complete the script
            if ((isFinal === true || String(isFinal) === 'true') && !currentSidebarItem) {
                currentSidebarItem = addCommandToSidebar('simplify');
            }
        }
    }

    console.log('=== 更新完成 ===\n');
}

// 解析并更新思考内容（保留作为备用）
function updateThinkingContent(content) {
    if (!content) return;

    console.log('收到思考内容:', content);

    // 提取各个字段
    const thinking = extractField(content, 'thinking');
    const command = extractField(content, 'command');
    const explanation = extractField(content, 'explanation');
    const isFinal = extractField(content, 'is_final');

    console.log('提取字段:', { thinking, command, explanation, isFinal });

    // 判断是否需要创建新卡片
    let needNewCard = false;

    if (!currentCard) {
        // 第一次，创建卡片
        needNewCard = true;
    } else if (thinking && lastThinkingValue) {
        // 检查是否是新一轮思考
        // 如果新内容比旧内容短很多，或者完全不同，说明是新一轮
        if (thinking.length < lastThinkingValue.length * 0.5) {
            needNewCard = true;
        } else if (!thinking.startsWith(lastThinkingValue.substring(0, 20)) &&
            !lastThinkingValue.startsWith(thinking.substring(0, 20))) {
            needNewCard = true;
        }
    }

    if (needNewCard) {
        console.log('创建新卡片');
        createNewCard();
    }

    // 在当前卡片内更新字段（流式更新）
    if (thinking) {
        updateField('thinking', thinking);
        lastThinkingValue = thinking;
        hasContent = true;
    }
    if (command) {
        updateField('command', command);
        hasContent = true;
    }
    if (explanation) {
        updateField('explanation', explanation);
        hasContent = true;
    }
    if (isFinal) {
        updateField('is_final', isFinal);
        hasContent = true;
    }
}

// 提取 JSON 字段
function extractField(content, fieldName) {
    try {
        // 先尝试解析完整JSON
        try {
            const jsonMatch = content.match(/\{[\s\S]*\}/);
            if (jsonMatch) {
                const json = JSON.parse(jsonMatch[0]);
                if (json[fieldName] !== undefined) {
                    return String(json[fieldName]);
                }
            }
        } catch (e) {
            // JSON解析失败，使用正则
        }

        // 使用正则提取（支持多行）
        // 模式1: "field": "value" (可能包含换行)
        let regex = new RegExp(`"${fieldName}":\\s*"([^"]*(?:\\\\.[^"]*)*)"`, 's');
        let match = content.match(regex);
        if (match) {
            return match[1].replace(/\\n/g, '\n').replace(/\\"/g, '"').trim();
        }

        // 模式2: "field": value (无引号，用于 boolean/number)
        regex = new RegExp(`"${fieldName}":\\s*([^,}\\s]+)`);
        match = content.match(regex);
        if (match) {
            return match[1].trim();
        }
    } catch (e) {
        console.error('提取字段失败:', fieldName, e);
    }
    return '';
}

// HTML 转义
function escapeHtml(text) {
    const div = document.createElement('div');
    div.textContent = text;
    return div.innerHTML;
}

// 示例卡片点击
const exampleCards = document.querySelectorAll('.example-card');
exampleCards.forEach(card => {
    card.addEventListener('click', () => {
        const expr = card.dataset.expr;
        const cond = card.dataset.cond;
        expressionInput.value = expr;
        conditionsInput.value = cond;
        expressionInput.focus();
    });
});

// 清空按钮
clearBtn.addEventListener('click', () => {
    expressionInput.value = '';
    conditionsInput.value = '';
    if (instructionInput) instructionInput.value = '';
    thinkingStepCount = 0;
    currentCard = null;
    currentSidebarItem = null;
    currentFields = {};
    lastThinkingValue = '';
    hasContent = false;
    lastStep = -1; // 重置步骤号

    // 重置所有滚动状态
    Object.values(scrollManagers).forEach(manager => manager.reset());

    resultOutput.innerHTML = '<div class="text-gray-400 text-center py-8">等待输入表达式...</div>';
    thinkingOutput.innerHTML = '<div class="text-gray-400 text-center py-8 text-xs">AI 思考过程将在这里显示...</div>';
    commandOutput.innerHTML = '<div class="text-gray-400 text-center py-8">执行的命令将在这里显示...</div>';
    errorOutput.innerHTML = '<div class="text-gray-400 text-center py-8">错误信息将在这里显示...</div>';
    stdCommandList.innerHTML = '<div class="text-gray-400 text-center py-8">求解命令将在此生成...</div>';
});

// 开始求解
solveBtn.addEventListener('click', async () => {
    console.log('求解按钮被点击');
    const expression = expressionInput.value.trim();
    console.log('表达式:', expression);

    if (!expression) {
        alert('请输入数学表达式');
        return;
    }

    if (isSolving) return;

    isSolving = true;
    updateSolveButton(true);

    // 重置思考步骤
    thinkingStepCount = 0;
    currentCard = null;
    currentSidebarItem = null; // Track current command item in sidebar
    currentFields = {};
    lastThinkingValue = '';
    hasContent = false;
    lastStep = -1; // 重置步骤号
    loadedSkillsThisSession.clear(); // Search-o1 风格：重置技能加载通知记录

    // 重置所有滚动状态
    Object.values(scrollManagers).forEach(manager => manager.reset());

    // 清空输出
    resultOutput.textContent = '';
    thinkingOutput.innerHTML = '';
    errorOutput.textContent = '';
    stdCommandList.innerHTML = ''; // Clear sidebar commands

    const conditions = conditionsInput.value.trim();

    // Intelligent initial command generation
    let initCmd = expression;
    const commandKeywords = ['calculate', 'prove', 'subgoal', 'defn', 'derive'];
    const firstWord = expression.split(' ')[0].toLowerCase();

    if (!commandKeywords.includes(firstWord)) {
        initCmd = `calculate ${expression}`;
    }

    if (conditions) {
        initCmd += ` for ${conditions}`;
    }
    addCommandToSidebar(initCmd);

    // 添加到历史记录
    addHistory(expression, conditions);

    try {
        // 使用 EventSource 进行流式输出
        const params = new URLSearchParams({
            expression: expression,
            conditions: conditions,
            instruction: instructionInput ? instructionInput.value.trim() : ''
        });

        eventSource = new EventSource(`/api/solve?${params.toString()}`);

        eventSource.onmessage = (event) => {
            const data = JSON.parse(event.data);

            console.log('收到数据:', data);

            // 更新各个输出区域
            if (data.commands !== undefined) {
                commandOutput.textContent = data.commands || '暂无命令';
                smartScroll(commandOutput);
            }
            if (data.results !== undefined) {
                resultOutput.textContent = data.results || '计算中...';
                smartScroll(resultOutput);
            }
            if (data.errors !== undefined) {
                if (data.errors && data.errors.trim().length > 0) {
                    errorOutput.textContent = data.errors;

                    // Remove current command if error occurs
                    if (currentSidebarItem) {
                        currentSidebarItem.remove();
                        currentSidebarItem = null;
                        console.log('Removed invalid command due to error');
                    }
                } else {
                    errorOutput.textContent = '暂无错误';
                }
                smartScroll(errorOutput);
            }
            if (data.thinking !== undefined && data.thinking) {
                console.log('思考内容:', data.thinking);
                // 现在thinking是一个对象，包含thinking, command, explanation, is_final
                updateThinkingFromObject(data.thinking);

                // Search-o1 风格：检测技能加载标记
                if (data.thinking.thinking) {
                    detectSkillLoading(data.thinking.thinking);
                }
            }
        };

        eventSource.onerror = (error) => {
            console.error('EventSource error:', error);
            eventSource.close();
            isSolving = false;
            updateSolveButton(false);

            if (resultOutput.textContent === '') {
                resultOutput.innerHTML = '<span class="text-red-500">连接失败，请检查后端服务是否运行</span>';
            }
        };

        eventSource.addEventListener('complete', () => {
            // 求解结束时，强制确标最后一步显示“证明结束”
            if (currentCard) {
                updateField('is_final', 'true');

                // Ensure script ends with a command (if last step didn't have one)
                if (!currentSidebarItem) {
                    addCommandToSidebar('simplify');
                }
            }

            eventSource.close();
            isSolving = false;
            updateSolveButton(false);
        });

    } catch (error) {
        console.error('Solve error:', error);
        errorOutput.textContent = `错误: ${error.message}`;
        isSolving = false;
        updateSolveButton(false);
    }
});

// 终止求解
stopBtn.addEventListener('click', async () => {
    if (!isSolving) return;

    try {
        await fetch('/api/stop', { method: 'POST' });

        if (eventSource) {
            eventSource.close();
            eventSource = null;
        }

        isSolving = false;
        updateSolveButton(false);

        resultOutput.textContent += '\n\n⚠ 已终止求解';
    } catch (error) {
        console.error('Stop error:', error);
    }
});

// 更新求解按钮状态
function updateSolveButton(solving) {
    if (solving) {
        solveBtnText.textContent = '⏳ 求解中...';
        solveSpinner.classList.remove('hidden');
        solveBtn.disabled = true;
        stopBtn.disabled = false;
    } else {
        solveBtnText.textContent = '🚀 开始求解';
        solveSpinner.classList.add('hidden');
        solveBtn.disabled = false;
        stopBtn.disabled = true;
    }
}

// 键盘快捷键

// Helper: Add command to sidebar
function addCommandToSidebar(text) {
    const div = document.createElement('div');
    div.className = 'group relative bg-gray-50 hover:bg-indigo-50 border border-gray-200 rounded p-2 cursor-pointer transition-all';
    div.onclick = function () { copyCommand(this, text); };

    div.innerHTML = `
        <code class="text-xs text-indigo-700 font-mono font-bold block">${text}</code>
        <span class="copy-hint absolute right-2 top-1.5 text-[10px] text-gray-400 opacity-0 group-hover:opacity-100 transition-opacity">点击复制</span>
    `;

    stdCommandList.appendChild(div);
    // Auto scroll list
    stdCommandList.scrollTop = stdCommandList.scrollHeight;

    return div;
}

// Helper: Update existing sidebar item
function updateCommandSidebarItem(element, text) {
    const code = element.querySelector('code');
    if (code) code.textContent = text;
    // Update onclick handler with new text
    element.onclick = function () { copyCommand(this, text); };
}


// 键盘快捷键
document.addEventListener('keydown', (e) => {
    // Ctrl/Cmd + Enter 开始求解
    if ((e.ctrlKey || e.metaKey) && e.key === 'Enter') {
        e.preventDefault();
        if (!isSolving) {
            solveBtn.click();
        }
    }

    // Esc 终止求解
    if (e.key === 'Escape' && isSolving) {
        stopBtn.click();
    }
});

// 复制所有命令
function copyAllCommands() {
    const commands = [];
    // Select all code blocks in the list
    const items = stdCommandList.querySelectorAll('code');
    items.forEach(code => {
        const text = code.innerText.trim();
        // Ignore placeholder
        if (text && !text.includes('求解命令将在此生成')) {
            commands.push(text);
        }
    });

    if (commands.length === 0) {
        // Visual feedback failure
        const btn = document.getElementById('copy-all-btn');
        btn.innerText = '❌ 无命令';
        setTimeout(() => btn.innerText = '📋 复制全部', 1000);
        return;
    }

    const text = commands.join('\n');
    navigator.clipboard.writeText(text).then(() => {
        const btn = document.getElementById('copy-all-btn');
        const originalText = "📋 复制全部";
        btn.innerText = '✅ 已复制';
        btn.classList.remove('text-indigo-600', 'bg-indigo-50');
        btn.classList.add('text-green-600', 'bg-green-50');

        setTimeout(() => {
            btn.innerText = originalText;
            btn.classList.add('text-indigo-600', 'bg-indigo-50');
            btn.classList.remove('text-green-600', 'bg-green-50');
        }, 1500);
    }).catch(err => {
        console.error('Copy failed', err);
    });
}

// 页面加载完成后获取当前模型
window.addEventListener('load', async () => {
    try {
        const response = await fetch('/api/config');
        const data = await response.json();
        if (data.model) {
            const modelEl = document.getElementById('current-model');
            if (modelEl) {
                modelEl.textContent = data.model;
            }
        }
        console.log('当前模型:', data.model);
    } catch (error) {
        console.error('Failed to fetch config:', error);
    }
});
