// 批量测试平台 - 前端逻辑

// ==================== 状态管理 ====================
const testState = {
    theories: [],           // 所有理论列表
    selectedTheories: [],   // 用户选择的理论
    isRunning: false,       // 是否正在测试
    summary: {              // 统计数据
        total: 0,
        success: 0,
        failed: 0,
        avgRounds: 0,
        totalTime: 0
    },
    problems: [],           // 所有题目的测试结果
    currentProblem: null,   // 当前正在测试的题目
    // Token 统计
    tokenStats: {
        totalCost: 0,
        totalTokens: 0,
        promptCacheHit: 0,
        promptCacheMiss: 0,
        completion: 0
    }
};

// ==================== DOM 元素 ====================
const elements = {
    // 左侧
    theoryList: document.getElementById('theoryList'),
    selectAllBtn: document.getElementById('selectAllBtn'),
    clearSelectionBtn: document.getElementById('clearSelectionBtn'),
    startTestBtn: document.getElementById('startTestBtn'),
    stopTestBtn: document.getElementById('stopTestBtn'),
    
    // 统计面板
    totalCount: document.getElementById('totalCount'),
    successCount: document.getElementById('successCount'),
    failedCount: document.getElementById('failedCount'),
    avgRounds: document.getElementById('avgRounds'),
    progressText: document.getElementById('progressText'),
    progressBar: document.getElementById('progressBar'),
    exportReportBtn: document.getElementById('exportReportBtn'),
    
    // 当前测试进度 - 已移除以避免卡死
    // currentTestPanel: document.getElementById('currentTestPanel'),
    // currentProblemId: document.getElementById('currentProblemId'),
    // currentRound: document.getElementById('currentRound'),
    // currentThinking: document.getElementById('currentThinking'),
    
    // 题目详情
    problemList: document.getElementById('problemList')
};

// ==================== 初始化 ====================
async function init() {
    await loadTheories();
    setupEventListeners();
}

// ==================== 加载理论列表 ====================
async function loadTheories() {
    try {
        const response = await fetch('/api/test/theories');
        if (!response.ok) {
            throw new Error('加载理论列表失败');
        }
        
        testState.theories = await response.json();
        renderTheoryList();
    } catch (error) {
        console.error('加载理论列表失败:', error);
        elements.theoryList.innerHTML = `
            <div class="text-center text-red-500 py-8">
                加载失败: ${error.message}
            </div>
        `;
    }
}

// ==================== 渲染理论列表 ====================
function renderTheoryList() {
    if (testState.theories.length === 0) {
        elements.theoryList.innerHTML = `
            <div class="text-center text-gray-400 py-8">
                未找到理论文件
            </div>
        `;
        return;
    }
    
    elements.theoryList.innerHTML = testState.theories.map(theory => `
        <label class="flex items-center gap-2 p-2 hover:bg-gray-50 rounded cursor-pointer">
            <input type="checkbox" 
                   class="theory-checkbox" 
                   value="${theory}"
                   ${testState.selectedTheories.includes(theory) ? 'checked' : ''}>
            <span class="text-sm text-gray-700">${theory}</span>
        </label>
    `).join('');
    
    // 添加复选框事件监听
    elements.theoryList.querySelectorAll('.theory-checkbox').forEach(checkbox => {
        checkbox.addEventListener('change', handleTheorySelection);
    });
}

// ==================== 处理理论选择 ====================
function handleTheorySelection(event) {
    const theory = event.target.value;
    if (event.target.checked) {
        if (!testState.selectedTheories.includes(theory)) {
            testState.selectedTheories.push(theory);
        }
    } else {
        testState.selectedTheories = testState.selectedTheories.filter(t => t !== theory);
    }
    updateStartButtonState();
}

// ==================== 更新开始按钮状态 ====================
function updateStartButtonState() {
    elements.startTestBtn.disabled = testState.selectedTheories.length === 0 || testState.isRunning;
}

// ==================== 事件监听器设置 ====================
function setupEventListeners() {
    // 全选
    elements.selectAllBtn.addEventListener('click', () => {
        testState.selectedTheories = [...testState.theories];
        renderTheoryList();
        updateStartButtonState();
    });
    
    // 清空选择
    elements.clearSelectionBtn.addEventListener('click', () => {
        testState.selectedTheories = [];
        renderTheoryList();
        updateStartButtonState();
    });
    
    // 开始测试
    elements.startTestBtn.addEventListener('click', startTest);
    
    // 停止测试
    elements.stopTestBtn.addEventListener('click', stopTest);
    
    // 导出报告
    elements.exportReportBtn.addEventListener('click', exportReport);
}

// ==================== 开始测试 ====================
function startTest() {
    if (testState.selectedTheories.length === 0) {
        alert('请至少选择一个理论文件');
        return;
    }
    
    // 重置状态
    testState.isRunning = true;
    testState.summary = { total: 0, success: 0, failed: 0, avgRounds: 0, totalTime: 0 };
    testState.problems = [];
    testState.currentProblem = null;
    
    // 更新 UI
    elements.startTestBtn.classList.add('hidden');
    elements.stopTestBtn.classList.remove('hidden');
    // elements.currentTestPanel.classList.remove('hidden');
    elements.exportReportBtn.classList.add('hidden');
    
    // 清空题目列表
    elements.problemList.innerHTML = '<div class="text-center text-gray-400 py-8">准备测试...</div>';
    
    // 启动 SSE 连接
    connectSSE();
}

// ==================== 停止测试 ====================
function stopTest() {
    testState.isRunning = false;
    elements.startTestBtn.classList.remove('hidden');
    elements.stopTestBtn.classList.add('hidden');
    
    // 关闭 SSE 连接
    if (testState.eventSource) {
        testState.eventSource.close();
        testState.eventSource = null;
    }
    
    // 发送停止请求
    fetch('/api/test/stop', { method: 'POST' })
        .catch(err => console.error('停止测试失败:', err));
}

// ==================== 更新 Token 显示 ====================
function updateTokenDisplay() {
    const totalCostEl = document.getElementById('totalCost');
    const totalTokensEl = document.getElementById('totalTokens');
    
    if (totalCostEl) {
        totalCostEl.textContent = testState.tokenStats.totalCost.toFixed(6);
    }
    if (totalTokensEl) {
        totalTokensEl.textContent = testState.tokenStats.totalTokens.toLocaleString();
    }
}

// ==================== 更新总结统计 ====================
function updateSummary() {
    const { total, success, failed } = testState.summary;
    const successRate = total > 0 ? ((success / total) * 100).toFixed(1) : 0;
    const progress = total > 0 ? ((success + failed) / total) * 100 : 0;
    
    // 计算平均轮数
    const totalRounds = testState.problems.reduce((sum, p) => sum + (p.num_rounds || 0), 0);
    const avgRounds = testState.problems.length > 0 ? (totalRounds / testState.problems.length).toFixed(1) : 0;
    
    elements.totalCount.textContent = total;
    elements.successCount.textContent = `${success} (${successRate}%)`;
    elements.failedCount.textContent = failed;
    elements.avgRounds.textContent = avgRounds;
    elements.progressText.textContent = `${progress.toFixed(0)}%`;
    elements.progressBar.style.width = `${progress}%`;
}

// ==================== 更新当前测试进度 ====================
// 已移除以避免卡死
// function updateCurrentTest(problemId, round, maxRound, thinking) {
//     testState.currentProblem = { problemId, round, maxRound, thinking };
//     
//     elements.currentProblemId.textContent = problemId || '-';
//     elements.currentRound.textContent = `${round}/${maxRound}`;
//     elements.currentThinking.textContent = thinking || '等待开始...';
// }

// ==================== 添加题目到列表 ====================
function addProblemToList(problem) {
    testState.problems.push(problem);
    renderProblemList();
    updateSummary();
}

// ==================== 更新题目状态 ====================
function updateProblemStatus(problemId, updates) {
    const problem = testState.problems.find(p => p.problem_id === problemId);
    if (problem) {
        Object.assign(problem, updates);
        renderProblemList();
        updateSummary();
    }
}

// ==================== 渲染题目列表 ====================
function renderProblemList() {
    if (testState.problems.length === 0) {
        elements.problemList.innerHTML = '<div class="text-center text-gray-400 py-8">暂无测试数据</div>';
        return;
    }
    
    elements.problemList.innerHTML = `
        <div class="space-y-2">
            ${testState.problems.map(problem => {
                const statusIcon = getStatusIcon(problem.status);
                const statusColor = getStatusColor(problem.status);
                return `
                    <div class="flex items-center justify-between p-3 bg-gray-50 rounded-lg hover:bg-gray-100 transition">
                        <div class="flex items-center gap-3 flex-1">
                            <span class="status-icon ${statusColor}">${statusIcon}</span>
                            <span class="font-mono text-sm text-gray-700">${problem.problem_id}</span>
                        </div>
                        <div class="flex items-center gap-4 text-sm text-gray-600">
                            <span>${problem.num_rounds || 0}轮</span>
                            <span>${problem.time_elapsed ? problem.time_elapsed.toFixed(1) + 's' : '-'}</span>
                        </div>
                    </div>
                `;
            }).join('')}
        </div>
    `;
}

// ==================== 获取状态图标 ====================
function getStatusIcon(status) {
    switch (status) {
        case 'success': return '✓';
        case 'failed': return '✗';
        case 'running': return '🔄';
        case 'waiting': return '⏳';
        default: return '⏳';
    }
}

// ==================== 获取状态颜色 ====================
function getStatusColor(status) {
    switch (status) {
        case 'success': return 'text-green-600';
        case 'failed': return 'text-red-600';
        case 'running': return 'text-indigo-600';
        case 'waiting': return 'text-gray-400';
        default: return 'text-gray-400';
    }
}

// ==================== 导出报告 ====================
function exportReport() {
    const report = generateMarkdownReport();
    const blob = new Blob([report], { type: 'text/markdown;charset=utf-8' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `test_report_${getTimestamp()}.md`;
    a.click();
    URL.revokeObjectURL(url);
}

// ==================== 生成 Markdown 报告 ====================
function generateMarkdownReport() {
    const { total, success, failed } = testState.summary;
    const successRate = total > 0 ? ((success / total) * 100).toFixed(1) : 0;
    const totalRounds = testState.problems.reduce((sum, p) => sum + (p.num_rounds || 0), 0);
    const avgRounds = testState.problems.length > 0 ? (totalRounds / testState.problems.length).toFixed(1) : 0;
    
    const successProblems = testState.problems.filter(p => p.status === 'success');
    const failedProblems = testState.problems.filter(p => p.status === 'failed');
    
    let report = `# 🧪 LLM-Iscalc 批量测试报告\n\n`;
    report += `**测试时间**: ${new Date().toLocaleString('zh-CN')}\n`;
    report += `**测试理论**: ${testState.selectedTheories.join(', ')}\n\n`;
    report += `---\n\n`;
    
    report += `## 📊 测试统计\n\n`;
    report += `| 指标 | 数值 |\n`;
    report += `|------|------|\n`;
    report += `| 总题目数 | ${total} |\n`;
    report += `| ✅ 成功 | ${success} (${successRate}%) |\n`;
    report += `| ❌ 失败 | ${failed} (${(100 - successRate).toFixed(1)}%) |\n`;
    report += `| 平均轮数 | ${avgRounds} |\n\n`;
    report += `---\n\n`;
    
    // Token & Cost 分析
    report += `## 💰 Token & Cost 分析\n\n`;
    report += `| 指标 | 总计 | 平均/题 |\n`;
    report += `|------|------|----------|\n`;
    report += `| **总成本** | **￥${testState.tokenStats.totalCost.toFixed(6)}** | **￥${(testState.tokenStats.totalCost / total).toFixed(6)}** |\n`;
    report += `| 总 Tokens | ${testState.tokenStats.totalTokens.toLocaleString()} | ${Math.round(testState.tokenStats.totalTokens / total).toLocaleString()} |\n`;
    report += `| Prompt (Hit) | ${testState.tokenStats.promptCacheHit.toLocaleString()} | ${Math.round(testState.tokenStats.promptCacheHit / total).toLocaleString()} |\n`;
    report += `| Prompt (Miss) | ${testState.tokenStats.promptCacheMiss.toLocaleString()} | ${Math.round(testState.tokenStats.promptCacheMiss / total).toLocaleString()} |\n`;
    report += `| Completion | ${testState.tokenStats.completion.toLocaleString()} | ${Math.round(testState.tokenStats.completion / total).toLocaleString()} |\n\n`;
    report += `---\n\n`;
    
    report += `## 📝 题目详情\n\n`;
    
    if (successProblems.length > 0) {
        report += `### ✅ 成功题目 (${successProblems.length}/${total})\n\n`;
        successProblems.forEach((p, idx) => {
            report += `#### ${idx + 1}. ${p.problem_id}\n`;
            report += `- **状态**: ✓ 成功\n`;
            report += `- **轮数**: ${p.num_rounds}\n`;
            report += `- **耗时**: ${p.time_elapsed ? p.time_elapsed.toFixed(1) + 's' : '-'}\n\n`;
            report += `---\n\n`;
        });
    }
    
    if (failedProblems.length > 0) {
        report += `### ❌ 失败题目 (${failedProblems.length}/${total})\n\n`;
        failedProblems.forEach((p, idx) => {
            report += `#### ${idx + 1}. ${p.problem_id}\n`;
            report += `- **状态**: ✗ 失败\n`;
            report += `- **轮数**: ${p.num_rounds}\n`;
            report += `- **耗时**: ${p.time_elapsed ? p.time_elapsed.toFixed(1) + 's' : '-'}\n\n`;
            report += `---\n\n`;
        });
    }
    
    report += `## 🔍 详细数据\n\n`;
    report += `| 题目ID | 状态 | 轮数 | 耗时(s) |\n`;
    report += `|--------|------|------|----------|\n`;
    testState.problems.forEach(p => {
        const status = p.status === 'success' ? '✓' : '✗';
        const time = p.time_elapsed ? p.time_elapsed.toFixed(1) : '-';
        report += `| ${p.problem_id} | ${status} | ${p.num_rounds} | ${time} |\n`;
    });
    
    report += `\n---\n\n`;
    report += `*报告生成时间: ${new Date().toLocaleString('zh-CN')}*\n`;
    
    return report;
}

// ==================== 获取时间戳 ====================
function getTimestamp() {
    const now = new Date();
    const year = now.getFullYear();
    const month = String(now.getMonth() + 1).padStart(2, '0');
    const day = String(now.getDate()).padStart(2, '0');
    const hour = String(now.getHours()).padStart(2, '0');
    const minute = String(now.getMinutes()).padStart(2, '0');
    const second = String(now.getSeconds()).padStart(2, '0');
    return `${year}${month}${day}_${hour}${minute}${second}`;
}

// ==================== 模拟测试（临时） ====================
function simulateTest() {
    // 模拟数据
    testState.summary.total = 15;
    
    const mockProblems = [
        { problem_id: 'standard_1', status: 'waiting', num_rounds: 0, time_elapsed: 0 },
        { problem_id: 'standard_2', status: 'waiting', num_rounds: 0, time_elapsed: 0 },
        { problem_id: 'standard_3', status: 'waiting', num_rounds: 0, time_elapsed: 0 },
        { problem_id: 'standard_4', status: 'waiting', num_rounds: 0, time_elapsed: 0 },
        { problem_id: 'standard_5', status: 'waiting', num_rounds: 0, time_elapsed: 0 }
    ];
    
    testState.problems = mockProblems;
    renderProblemList();
    updateSummary();
    
    // 模拟测试过程
    let currentIdx = 0;
    const interval = setInterval(() => {
        if (!testState.isRunning || currentIdx >= mockProblems.length) {
            clearInterval(interval);
            // elements.currentTestPanel.classList.add('hidden');
            elements.exportReportBtn.classList.remove('hidden');
            elements.startTestBtn.classList.remove('hidden');
            elements.stopTestBtn.classList.add('hidden');
            return;
        }
        
        const problem = mockProblems[currentIdx];
        const isSuccess = Math.random() > 0.3;
        const numRounds = Math.floor(Math.random() * 8) + 2;
        const timeElapsed = Math.random() * 5 + 0.5;
        
        // 更新当前测试
        // updateCurrentTest(problem.problem_id, numRounds, 10, `正在求解 ${problem.problem_id}...\napply integral identity\nsimplify`);
        
        // 更新题目状态
        updateProblemStatus(problem.problem_id, {
            status: 'running',
            num_rounds: numRounds
        });
        
        // 1秒后完成
        setTimeout(() => {
            updateProblemStatus(problem.problem_id, {
                status: isSuccess ? 'success' : 'failed',
                num_rounds: numRounds,
                time_elapsed: timeElapsed
            });
            
            if (isSuccess) {
                testState.summary.success++;
            } else {
                testState.summary.failed++;
            }
            
            updateSummary();
            currentIdx++;
        }, 1000);
        
    }, 2000);
}

// ==================== SSE 连接 ====================
function connectSSE() {
    const theories = testState.selectedTheories.join(',');
    const url = `/api/test/run?theories=${encodeURIComponent(theories)}&max_workers=5&max_step=10`;
    
    testState.eventSource = new EventSource(url);
    
    // 测试开始
    testState.eventSource.addEventListener('test_start', (event) => {
        const data = JSON.parse(event.data);
        testState.summary.total = data.total_problems;
        updateSummary();
        console.log('测试开始:', data);
    });
    
    // 题目开始
    testState.eventSource.addEventListener('problem_start', (event) => {
        const data = JSON.parse(event.data);
        const problem = {
            problem_id: data.problem_id,
            problem_text: data.problem_text,
            status: 'running',
            num_rounds: 0,
            time_elapsed: 0
        };
        addProblemToList(problem);
        // updateCurrentTest(data.problem_id, 0, 10, '准备开始...');
        console.log('题目开始:', data);
    });
    
    // round_update 事件已移除以避免卡死
    
    // AI 思考流式输出 - 已移除以避免卡死
    // testState.eventSource.addEventListener('thinking_stream', (event) => {
    //     const data = JSON.parse(event.data);
    //     const currentThinking = elements.currentThinking.textContent;
    //     elements.currentThinking.textContent = currentThinking + data.chunk;
    //     // 自动滚动到底部
    //     elements.currentThinking.scrollTop = elements.currentThinking.scrollHeight;
    // });
    
    // 题目完成
    testState.eventSource.addEventListener('problem_complete', (event) => {
        const data = JSON.parse(event.data);
        updateProblemStatus(data.problem_id, {
            status: data.status,
            num_rounds: data.num_rounds,
            time_elapsed: data.time_elapsed
        });
        
        // 更新统计
        if (data.status === 'success') {
            testState.summary.success++;
        } else {
            testState.summary.failed++;
        }
        updateSummary();
        
        // 累加 Token 统计
        if (data.token_usage) {
            testState.tokenStats.totalCost += data.token_usage.cost;
            testState.tokenStats.totalTokens += data.token_usage.total_tokens;
            testState.tokenStats.promptCacheHit += data.token_usage.prompt_cache_hit;
            testState.tokenStats.promptCacheMiss += data.token_usage.prompt_cache_miss;
            testState.tokenStats.completion += data.token_usage.completion;
            
            // 更新显示
            updateTokenDisplay();
        }
        
        // 清空当前测试面板
        // updateCurrentTest(null, 0, 10, '等待下一题...');
        console.log('题目完成:', data);
    });
    
    // 测试完成
    testState.eventSource.addEventListener('test_complete', (event) => {
        console.log('测试完成');
        testState.isRunning = false;
        testState.eventSource.close();
        testState.eventSource = null;
        
        // 隐藏当前测试面板，显示导出按钮
        // elements.currentTestPanel.classList.add('hidden');
        elements.exportReportBtn.classList.remove('hidden');
        elements.startTestBtn.classList.remove('hidden');
        elements.stopTestBtn.classList.add('hidden');
    });
    
    // 错误处理
    testState.eventSource.addEventListener('error', (event) => {
        console.error('SSE 错误:', event);
        if (testState.eventSource.readyState === EventSource.CLOSED) {
            console.log('SSE 连接已关闭');
            testState.isRunning = false;
            elements.startTestBtn.classList.remove('hidden');
            elements.stopTestBtn.classList.add('hidden');
        }
    });
    
    testState.eventSource.onerror = (error) => {
        console.error('SSE 连接错误:', error);
    };
}

// ==================== 页面加载完成后初始化 ====================
document.addEventListener('DOMContentLoaded', init);
