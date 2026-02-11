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
// const toggleMaximizeBtn = document.getElementById('toggle-maximize-btn'); // Removed specific ID
const thinkingSection = document.getElementById('thinking-section');
const historyList = document.getElementById('history-list');
const clearHistoryBtn = document.getElementById('clear-history-btn');
const maximizeBackdrop = document.getElementById('maximize-backdrop');
let currentMaximizedSection = null;

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

function createHistoryItemDOM(item) {
    const div = document.createElement('div');
    div.className = 'history-item';
    // 设置唯一 ID (用于动画追踪)
    div.dataset.key = `${item.expression}||${item.conditions || ''}`;

    div.innerHTML = `
        <div class="history-item-expr" title="${item.expression}">${item.expression}</div>
        ${item.conditions ? `<div class="history-item-cond" title="${item.conditions}">条件: ${item.conditions}</div>` : ''}
        <div class="history-item-time">${item.timestamp}</div>
        <button class="history-delete-btn" title="删除">🗑️</button>
    `;

    div.addEventListener('click', () => {
        expressionInput.value = item.expression;
        conditionsInput.value = item.conditions;
        expressionInput.focus();
    });

    const deleteBtn = div.querySelector('.history-delete-btn');
    if (deleteBtn) {
        deleteBtn.addEventListener('click', (e) => {
            e.stopPropagation();
            deleteHistoryItem(item);
        });
    }

    return div;
}

// 删除历史记录项
function deleteHistoryItem(itemToDelete) {
    if (!confirm('确定要删除这条记录吗？')) return;

    // Update data model
    const index = solveHistory.findIndex(h =>
        h.expression === itemToDelete.expression &&
        h.conditions === itemToDelete.conditions &&
        h.timestamp === itemToDelete.timestamp
    );

    if (index !== -1) {
        solveHistory.splice(index, 1);
        saveHistory();

        // Update DOM
        const key = `${itemToDelete.expression}||${itemToDelete.conditions || ''}`;

        const items = historyList.querySelectorAll('.history-item');
        let el = null;
        items.forEach(item => {
            if (item.dataset.key === key) el = item;
        });

        if (el) {
            // Animation for removal
            el.style.transform = 'translateX(100px)';
            el.style.opacity = '0';
            setTimeout(() => {
                if (solveHistory.length === 0) {
                    renderHistory(); // Show empty message
                } else {
                    el.remove();
                }
            }, 300);
        } else {
            renderHistory(); // Fallback
        }
    }
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

    /**
     * Setup maximize functionality with FLIP animation
     * @param {HTMLElement} btn - The toggle button
     * @param {HTMLElement} section - The section to maximize
     */
    function setupMaximize(btn, section) {
        btn.addEventListener('click', () => {
            // 1. First: Record start state
            const firstRect = section.getBoundingClientRect();
            const isMaximized = section.classList.contains('maximized');

            // 2. Last: Apply class to change layout
            section.classList.toggle('maximized');

            // 3. Invert: Measure new position and calculate delta
            const lastRect = section.getBoundingClientRect();
            const deltaX = firstRect.left - lastRect.left;
            const deltaY = firstRect.top - lastRect.top;
            const deltaW = firstRect.width / lastRect.width;
            const deltaH = firstRect.height / lastRect.height;

            // Apply transform to make it look like it's still at the start position
            section.style.transformOrigin = 'top left';
            section.style.transform = `translate(${deltaX}px, ${deltaY}px) scale(${deltaW}, ${deltaH})`;
            section.style.transition = 'none'; // Disable transition for instant setup

            // Trigger reflow
            section.offsetHeight;

            // 4. Play: Enable transition and remove transform to animate to end state
            section.style.transition = ''; // Restore CSS transition
            section.style.transform = '';

            // Update button icon
            btn.textContent = isMaximized ? '⤢' : '⤡';
        });
    }
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


// SVG Icons
const ICON_MAXIMIZE = `
<svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke-width="1.5" stroke="currentColor" class="w-4 h-4">
  <path stroke-linecap="round" stroke-linejoin="round" d="M3.75 3.75v4.5m0-4.5h4.5m-4.5 0L9 9M3.75 20.25v-4.5m0 4.5h4.5m-4.5 0L9 15M20.25 3.75h-4.5m4.5 0v4.5m0-4.5L15 9m5.25 11.25h-4.5m4.5 0v-4.5m0 4.5L15 15" />
</svg>`;

const ICON_RESTORE = `
<svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke-width="1.5" stroke="currentColor" class="w-4 h-4">
  <path stroke-linecap="round" stroke-linejoin="round" d="M9 9V4.5M9 9H4.5M9 9L3.75 3.75M9 15v4.5M9 15H4.5M9 15l-5.25 5.25M15 9h4.5M15 9V4.5M15 9l5.25-5.25M15 15h4.5M15 15v4.5M15 15l5.25 5.25" />
</svg>`;

// Initialize generic maximize buttons
document.querySelectorAll('.maximize-btn').forEach(btn => {
    // Find target section: parent with class ending in '-section'
    const section = btn.closest('.thinking-section, .std-command-section, .result-section, .command-section, .error-section');
    if (section) {
        setupMaximize(btn, section);
    }
});

/**
 * Setup maximize functionality with Portal Strategy and Layout Animation (No Distortion)
 * Moves the section to body when maximized and animates width/height/top/left
 * @param {HTMLElement} btn - The toggle button
 * @param {HTMLElement} section - The section to maximize
 */
function setupMaximize(btn, section) {
    // Store placeholder reference on the section to retrieve it later
    section._placeholder = null;

    btn.addEventListener('click', () => {
        const isMaximized = section.classList.contains('maximized');

        if (!isMaximized) {
            // === MAXIMIZE ACTION ===

            // 1. Record start state (block layout)
            const firstRect = section.getBoundingClientRect();

            // 2. Create Placeholder
            const placeholder = document.createElement('div');
            placeholder.style.width = firstRect.width + 'px';
            placeholder.style.height = firstRect.height + 'px';
            placeholder.style.flexGrow = window.getComputedStyle(section).flexGrow;
            placeholder.style.flexShrink = window.getComputedStyle(section).flexShrink;
            placeholder.className = section.className.replace('maximized', '');
            placeholder.style.visibility = 'hidden';
            placeholder.id = 'placeholder-' + Math.random().toString(36).substr(2, 9);

            // 3. Portal Move
            section.parentNode.insertBefore(placeholder, section);
            section._placeholder = placeholder;
            document.body.appendChild(section);

            // 4. Set Start State (Fixed at original position)
            section.classList.add('maximized'); // Add class for base styles (e.g. z-index if needed, though we set manually)
            section.style.position = 'fixed';
            section.style.top = firstRect.top + 'px';
            section.style.left = firstRect.left + 'px';
            section.style.width = firstRect.width + 'px';
            section.style.height = firstRect.height + 'px';
            section.style.margin = '0';
            section.style.zIndex = '1000';
            section.style.transform = '';
            section.style.transition = 'none'; // Instant placement

            // Lock Content Width to prevent reflow
            // Only lock for .output-box (Result/Command/Error) to prevent text reflow.
            // Thinking Section (#thinking-output) is Intentionally Excluded so it expands to full width.
            const content = section.querySelector('.output-box');
            if (content) {
                // Lock to current offsetWidth to maintain text wrap state
                content.style.width = content.offsetWidth + 'px';
                // Also ensure flex-grow doesn't stretch it immediately if we change container
                content.style.flex = 'none';
                // Align to Top-Left (User Request)
                content.style.margin = '0';
            }

            // Force Reflow
            section.offsetHeight;

            // 5. Animate to End State (Fullscreen-ish)
            // Use Layout Animation instead of Transform to avoid text distortion
            section.style.transition = 'all 0.3s cubic-bezier(0.4, 0, 0.2, 1)';
            section.style.top = '20px';
            section.style.left = '20px';
            section.style.width = 'calc(100vw - 40px)';
            section.style.height = 'calc(100vh - 40px)';

            // Unlock Content Width after animation to fill space
            const onMaximizeEnd = () => {
                if (content) {
                    content.style.width = '';
                    content.style.flex = ''; // Revert to stylesheet default (flex: 1)
                }
                section.removeEventListener('transitionend', onMaximizeEnd);
            };
            section.addEventListener('transitionend', onMaximizeEnd);

            // Update UI
            btn.innerHTML = ICON_RESTORE;
            maximizeBackdrop.classList.remove('pointer-events-none', 'opacity-0');
            maximizeBackdrop.classList.add('pointer-events-auto', 'opacity-100');
            currentMaximizedSection = { btn, section };
            toggleScrollLock(true);

        } else {
            // === RESTORE ACTION ===

            const placeholder = section._placeholder;
            if (!placeholder) return;

            // 1. Get current Fixed state
            const rect = section.getBoundingClientRect();
            // 2. Get target state (Placeholder)
            const targetRect = placeholder.getBoundingClientRect();

            // 3. Set Start positions explicitly (to ensure transition starts from current)
            section.style.top = rect.top + 'px';
            section.style.left = rect.left + 'px';
            section.style.width = rect.width + 'px';
            section.style.height = rect.height + 'px';
            // Ensure transition is active
            section.style.transition = 'all 0.3s cubic-bezier(0.4, 0, 0.2, 1)';

            section.offsetHeight; // Force reflow

            // 4. Animate to Target
            section.style.top = targetRect.top + 'px';
            section.style.left = targetRect.left + 'px';
            section.style.width = targetRect.width + 'px';
            section.style.height = targetRect.height + 'px';

            // UI Updates
            btn.innerHTML = ICON_MAXIMIZE;
            maximizeBackdrop.classList.add('pointer-events-none', 'opacity-0');
            maximizeBackdrop.classList.remove('pointer-events-auto', 'opacity-100');
            currentMaximizedSection = null;
            toggleScrollLock(false);

            // 5. Cleanup logic
            const onTransitionEnd = () => {
                if (section.style.width !== targetRect.width + 'px') {
                    // Safety check: ensure we don't trigger if interrupted (though unlikely with this logic)
                }

                section.classList.remove('maximized');
                section.style.cssText = ''; // Clear all inline styles

                // Unlock Content Width
                const content = section.querySelector('.output-box');
                if (content) {
                    content.style.width = '';
                    content.style.flex = '';
                    content.style.margin = '';
                }

                // Portal Back
                if (placeholder.parentNode) {
                    placeholder.parentNode.insertBefore(section, placeholder);
                    placeholder.remove();
                }
                section._placeholder = null;

                // Clean up listener
                section.removeEventListener('transitionend', onTransitionEnd);
            };
            section.addEventListener('transitionend', onTransitionEnd); // Only trigger once
        }
    }); // End click listener
}

/**
 * Toggle scroll lock on body to prevent background scrolling when modal/maximized view is active
 * Handles scrollbar width compensation to prevent layout shift
 * @param {boolean} enable - true to lock, false to unlock
 */
function toggleScrollLock(enable) {
    const body = document.body;
    if (enable) {
        // Calculate scrollbar width
        const scrollbarWidth = window.innerWidth - document.documentElement.clientWidth;
        // Add padding to prevent layout shift
        body.style.paddingRight = `${scrollbarWidth}px`;
        body.style.overflow = 'hidden';
    } else {
        body.style.paddingRight = '';
        body.style.overflow = '';
    }
}

// Backdrop click to minimize
if (maximizeBackdrop) {
    maximizeBackdrop.addEventListener('click', () => {
        if (currentMaximizedSection) {
            // Simulate click on the toggle button of the currently maximized section
            currentMaximizedSection.btn.click();
        }
    });
}


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

// 侧边栏同步状态
let lastSyncedResultLineCount = 0; // Track how many lines we've processed
let lastSyncedCommandCount = 0; // Track how many commands we've synced

/**
 * Sync sidebar commands with lines in the result output
 * Parses lines like: "= ... (command)" and adds them to sidebar if missing
 * @param {string} resultText - Full text of the result output
 */
function syncSidebarWithResults(resultText) {
    if (!resultText) return;

    // 2026-02-01 再次优化：改进正则以处理跨行和复杂命令格式
    // 策略：分两步
    // 1. 先找出所有以 = 开头的块（可能跨多行）
    // 2. 从每个块中提取最后一个括号内的内容

    // 按行分割，然后找出所有以 = 开头的行
    const lines = resultText.split('\n');
    const commandsFound = [];

    for (let i = 0; i < lines.length; i++) {
        const line = lines[i].trim();
        // 只处理以 = 开头的行
        if (line.startsWith('=')) {
            // 提取这一行中所有的 (xxx) 内容
            // 使用更智能的括号匹配：找到最后一个完整的括号对
            const bracketMatches = [];
            let depth = 0;
            let start = -1;

            for (let j = 0; j < line.length; j++) {
                if (line[j] === '(') {
                    if (depth === 0) start = j;
                    depth++;
                } else if (line[j] === ')') {
                    depth--;
                    if (depth === 0 && start !== -1) {
                        // 找到一个完整的括号对
                        const content = line.substring(start + 1, j).trim();
                        bracketMatches.push(content);
                        start = -1;
                    }
                }
            }

            // 取最后一个括号内容作为命令
            if (bracketMatches.length > 0) {
                const lastCommand = bracketMatches[bracketMatches.length - 1];
                // 清理命令：去掉可能的前缀标记（如 "? "）
                const cleanedCommand = lastCommand.replace(/^\?\s*/, '').trim();
                if (cleanedCommand) {
                    commandsFound.push(cleanedCommand);
                }
            }
        }
    }

    // 过滤掉无效命令和提示性内容
    // 用户可根据需要调整此列表
    const invalidCommands = [
        'calculate',     // 初始命令已在顶部显示
        'LLM判定',       // 提示性内容，非实际命令
        'LLM判断',       // 同上
        '计算',          // 同上
        '{calculate}',   // 格式化标记
        'simplify'       // 用户要求：不显示最后的 simplify（2026-02-01）
    ];
    commandsFound = commandsFound.filter(cmd => !invalidCommands.includes(cmd));

    // Iterate through new commands we haven't processed yet
    // We assume resultText is strictly append-only so we can skip the first N commands we already synced
    // However, to be safe against overlaps (if thinking added one), we check the last item.

    // Strategy:
    // We want the sidebar to reflect 'commandsFound'.
    // 'commandsFound' starts with the first step.
    // The sidebar has [Initial Command, Step 1, Step 2...]
    // We skip 'lastSyncedCommandCount' items from 'commandsFound' to only process new ones.

    const newCommands = commandsFound.slice(lastSyncedCommandCount);

    newCommands.forEach(cmd => {
        // Check if this command was just added by the "Thinking" stream
        // Get last item in sidebar
        const lastSidebarItem = stdCommandList.lastElementChild;
        let lastSidebarText = '';
        if (lastSidebarItem) {
            const code = lastSidebarItem.querySelector('code');
            if (code) lastSidebarText = code.innerText.trim();
        }

        // Deduplicate: If the last sidebar item is identical to this result command, assume it's the same step
        if (lastSidebarText === cmd) {
            // Already there (likely from thinking stream), just mark as synced
            console.log('Sidebar sync: Skipping duplicate', cmd);
        }
        // Also check if lastSidebarText is "simplify" and cmd is "simplify" (common case)
        else {
            // Not in sidebar (or different), so add it
            // This catches hidden intermediate steps like "partial-fraction"
            addCommandToSidebar(cmd);
            console.log('Sidebar sync: Added missing command', cmd);

            // If we added it, it becomes the new "currentSidebarItem" context potentially, 
            // but usually we don't need to link it to currentCard unless we want to update it later.
            // For now, standalone addition is fine.
        }
    });

    // Update counter
    lastSyncedCommandCount = commandsFound.length;
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

// 添加命令到侧边栏
function addCommandToSidebar(command, isThinking = false, step = 0) {
    if (!command) return;

    // Filter out "read_skill" commands from the sidebar as requested
    if (command.trim().startsWith('read_skill')) {
        return;
    }

    // 检查是否重复（最后一条如果一样就不加）
    const lastSidebarItem = stdCommandList.lastElementChild;
    let lastSidebarText = '';
    if (lastSidebarItem) {
        const code = lastSidebarItem.querySelector('code');
        if (code) lastSidebarText = code.innerText.trim();
    }

    // Skip duplicates unless it's a new step (to avoid merging different steps)
    // But mostly we want to avoid duplicates regardless
    if (lastSidebarText === command.trim()) {
        console.log('Sidebar: Skipping duplicate command:', command);
        return lastSidebarItem; // 返回已存在的元素
    }

    const commandItem = document.createElement('div');
    commandItem.className = `sidebar-command-item ${isThinking ? 'thinking-command' : ''}`;
    // Assign data-step attribute if step is valid
    if (step > 0) {
        commandItem.dataset.step = step;
    }
    commandItem.innerHTML = `<code>${command}</code>`;
    stdCommandList.appendChild(commandItem);
    smartScroll(stdCommandList); // 滚动侧边栏到底部
    return commandItem;
}

// 更新侧边栏命令项
function updateCommandSidebarItem(item, newCommand) {
    if (!item || !newCommand) return;

    // Filter out "read_skill" commands from the sidebar as requested
    if (newCommand.trim().startsWith('read_skill')) {
        // If the item was previously a valid command and now it's read_skill,
        // we might want to remove it or just not update it.
        // For now, let's just not update it.
        return;
    }

    const code = item.querySelector('code');
    if (code) {
        code.innerText = newCommand;
    }
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
            // 2026-02-01 重大修改：完全禁用从 thinking 向 sidebar 添加命令的逻辑
            // 原因：用户明确要求"命令框中的命令和结果框中的结果括号中的命令对应起来"
            // 解决方案：让 syncSidebarWithResults 成为 sidebar 的唯一数据源
            // 这样可以保证左右两边 100% 一致，不再出现乱序或不匹配的情况
            /*
            // Filter out "None" or empty commands
            if (command && command !== "None" && command.trim() !== "") {
                // 注意:不再区分用户命令和中间步骤,所有命令都添加
                // 因为中间步骤现在会通过独立的COMMAND事件发送
                if (!currentSidebarItem) {
                    currentSidebarItem = addCommandToSidebar(command, true, step); // Pass step!
                } else {
                    updateCommandSidebarItem(currentSidebarItem, command);
                }
            }
            */

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

            // 2026-02-01: 禁用此逻辑，避免与 syncSidebarWithResults 冲突
            // If final and no command yet for this step, add 'simplify' to complete the script
            // if ((isFinal === true || String(isFinal) === 'true') && !currentSidebarItem) {
            //     currentSidebarItem = addCommandToSidebar('simplify');
            // }
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

    // Reset sync state
    lastSyncedResultLineCount = 0;
    lastSyncedCommandCount = 0;

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
        // 检查是否为等式：包含等号 '='
        // 简单检测：如果表达式中包含 '=' 且不是以 '<' 或 '>' 结尾（排除 <= 和 >=）
        const hasEquality = expression.includes('=') &&
            !expression.includes('<=') &&
            !expression.includes('>=') &&
            !expression.includes('!=');

        if (hasEquality) {
            initCmd = `prove ${expression}`;
        } else {
            initCmd = `calculate ${expression}`;
        }
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

                // 2026-02-01 架构重构：移除正则匹配逻辑
                // 不再从results解析命令，直接监听 COMMAND_SUCCESS 事件
                // syncSidebarWithResults(data.results);  // 已删除
            }

            // 2026-02-01 架构重构：监听命令执行成功事件
            if (data.command_success) {
                const cmd = data.command_success;
                const step = cmd.step || 0;
                const commandText = cmd.content || '';

                // 过滤不需要显示的命令
                // 过滤不需要显示的命令
                const excludeCommands = ['read_skill'];
                // const excludeCommands = ['read_skill', 'simplify']; // 2026-02-01: 恢复 simplify 显示，否则中间步骤也不显示了
                const shouldAdd = !excludeCommands.some(exc => commandText.trim().startsWith(exc));

                if (shouldAdd && commandText) {
                    addCommandToSidebar(commandText, false, step);
                    console.log('✅ 命令执行成功，已添加:', commandText);
                }
            }

            // 2026-02-01 架构重构：监听命令执行失败事件
            if (data.command_failure) {
                const cmd = data.command_failure;
                const step = cmd.step || 0;
                const commandText = cmd.content || '';

                // 命令失败，不添加到 sidebar
                console.log('❌ 命令执行失败，不添加:', commandText);

                // 如果该命令已经被添加（例如在 thinking 阶段），删除它
                const existingItem = stdCommandList.querySelector(`.sidebar-command-item[data-step="${step}"]`);
                if (existingItem) {
                    existingItem.remove();
                    console.log('❌ 已删除失败命令:', commandText);
                }
            }

            if (data.errors !== undefined) {
                if (data.errors && data.errors.trim().length > 0) {
                    errorOutput.textContent = data.errors;

                    // 2026-02-01 架构重构：移除复杂的错误匹配删除逻辑
                    // 命令的添加/删除现在完全由 COMMAND_SUCCESS/COMMAND_FAILURE 事件控制
                    // 不再需要解析错误信息中的步骤号
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

        // Backend signals completion explicitly via SSE event: "complete"
        // Keep this handler minimal: close stream + unlock UI.
        eventSource.addEventListener('complete', () => {
            try {
                eventSource.close();
            } catch (e) {
                // ignore
            }
            eventSource = null;
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
