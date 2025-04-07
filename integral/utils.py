"""Utility functions"""

def indent(s: str, n: int = 2) -> str:
    """Indent input string by given number of spaces."""
    lines = s.split('\n')
    lines = [' ' * n + line for line in lines]
    return '\n'.join(lines)
