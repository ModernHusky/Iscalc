"""Session Manager 模块

管理求解会话和历史。
"""

import uuid
from dataclasses import dataclass, field
from datetime import datetime
from typing import Optional, Dict, List, Any


@dataclass
class Session:
    """会话数据结构"""
    id: str
    expression: str
    history: List[Dict[str, Any]] = field(default_factory=list)
    status: str = "active"  # active, completed, error
    created_at: datetime = field(default_factory=datetime.now)
    updated_at: datetime = field(default_factory=datetime.now)
    result: Optional[Any] = None
    final_expression: Optional[str] = None
    latex_solution: Optional[str] = None


class SessionManager:
    """会话管理器"""
    
    def __init__(self):
        self.sessions: Dict[str, Session] = {}
    
    def create_session(self, expression: str) -> Session:
        """创建新会话"""
        session_id = str(uuid.uuid4())[:8]
        session = Session(
            id=session_id,
            expression=expression
        )
        self.sessions[session_id] = session
        return session
    
    def get_session(self, session_id: str) -> Optional[Session]:
        """获取会话"""
        return self.sessions.get(session_id)
    
    def update_session(
        self,
        session_id: str,
        history: Optional[List[Dict[str, Any]]] = None,
        status: Optional[str] = None,
        result: Optional[Any] = None,
        final_expression: Optional[str] = None,
        latex_solution: Optional[str] = None
    ) -> Optional[Session]:
        """更新会话"""
        session = self.sessions.get(session_id)
        if session is None:
            return None
        
        if history is not None:
            session.history = history
        if status is not None:
            session.status = status
        if result is not None:
            session.result = result
        if final_expression is not None:
            session.final_expression = final_expression
        if latex_solution is not None:
            session.latex_solution = latex_solution
        
        session.updated_at = datetime.now()
        return session
    
    def add_step(self, session_id: str, step: Dict[str, Any]) -> Optional[Session]:
        """添加步骤到会话历史"""
        session = self.sessions.get(session_id)
        if session is None:
            return None
        
        session.history.append(step)
        session.updated_at = datetime.now()
        return session
    
    def reset_session(self, session_id: str) -> Optional[Session]:
        """重置会话"""
        session = self.sessions.get(session_id)
        if session is None:
            return None
        
        session.history = []
        session.status = "active"
        session.result = None
        session.final_expression = None
        session.latex_solution = None
        session.updated_at = datetime.now()
        return session
    
    def delete_session(self, session_id: str) -> bool:
        """删除会话"""
        if session_id in self.sessions:
            del self.sessions[session_id]
            return True
        return False
    
    def list_sessions(self) -> List[Session]:
        """列出所有会话"""
        return list(self.sessions.values())
    
    def get_active_sessions(self) -> List[Session]:
        """获取所有活跃会话"""
        return [s for s in self.sessions.values() if s.status == "active"]
