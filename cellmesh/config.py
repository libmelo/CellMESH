"""
CELL MESH 配置文件
集中管理所有参数、阈值和路径设置
"""
from pathlib import Path
from typing import Literal, Dict, Any

# ==================== 基础路径配置 ====================
# 数据目录
DATA_DIR = Path(__file__).parent / "data"

# 确保目录存在
DATA_DIR.mkdir(exist_ok=True)

# ==================== 通用阈值配置 ====================
# 最小细胞数阈值，低于该数量的细胞类型会被自动排除
MIN_CELL_COUNT: int = 100

# 最小细胞表达比例阈值
MIN_EXPR_FRAC: float = 0.05

# ==================== 代谢物可用性计算默认参数 ====================
METABOLITE_AVAILABILITY_DEFAULTS: Dict[str, Any] = {
    "lower": 5,
    "upper": 95,
    "eps": 0.05,
    "beta": 0.5,
    "missing_C_norm": 0.2,
    "missing_E_norm": 0.5,
    "min_cells": 1,
}

# ==================== 角色和传感器类型常量 ====================
# 合法的酶角色类型
VALID_ROLES = {"production", "degradation", "export"}

# 合法的传感器类型
VALID_SENSOR_TYPES = {"surface_receptor", "transporter", "nuclear_receptor", "intracellular_sensor"}

# 角色到反应方向的映射（新 availability 算法使用）
ROLE_TO_DIRECTION = {
    'production': 'product',      # 产生代谢物 → P 矩阵
    'degradation': 'substrate',    # 降解代谢物 → C 矩阵
    'export': 'exporter',          # 外排代谢物 → E 矩阵
}
