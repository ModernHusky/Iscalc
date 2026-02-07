
import sys
import os

# Ensure we can import from the project root
sys.path.append(os.getcwd())

from llm_iscalc.command_executor import CommandExecutor

def test_is_finished_bug():
    print("=== Testing is_finished bug ===")
    
    executor = CommandExecutor()
    
    # 1. Initialize
    expr = "prove (INT x:[1,oo]. 1/(x*(x^2+1))) = log(2)/2"
    print(f"\n[1] Initializing with: {expr}")
    res = executor.initialize(expr, [])
    if not res.success:
        print(f"FAILED to initialize: {res.error}")
        return
        
    state_name = executor.get_current_state_name()
    finished = executor.is_finished()
    print(f"State: {state_name}")
    print(f"is_finished: {finished}")
    
    if finished:
        print("❌ BUG: Finished immediately after initialization!")
    
    # 2. Execute lhs:
    print(f"\n[2] Executing 'lhs:'")
    res = executor.execute("lhs:")
    if not res.success:
        print(f"FAILED to execute lhs:: {res.error}")
        return
        
    state_name = executor.get_current_state_name()
    finished = executor.is_finished()
    print(f"State: {state_name}")
    print(f"is_finished: {finished}")
    
    if finished:
        print("❌ BUG: Finished after 'lhs:' but proof is NOT done!")
    else:
        print("✅ Correct: Not finished after 'lhs:'")

if __name__ == "__main__":
    test_is_finished_bug()
