"""
评分相关工具函数
"""
import numpy as np
import pandas as pd


def sigmoid(x: np.ndarray) -> np.ndarray:
    """
    Sigmoid 激活函数，将值映射到 (0, 1) 范围
    
    参数:
        x: 输入数组
    
    返回:
        Sigmoid 转换后的数组
    """
    return 1 / (1 + np.exp(-x))


def zscore_by_gene(expr_matrix: pd.DataFrame) -> pd.DataFrame:
    """
    对每个基因（列）进行 z-score 标准化
    
    参数:
        expr_matrix: 表达矩阵，行为样本/细胞类型，列为基因
    
    返回:
        z-score 标准化后的矩阵
    """
    mean = expr_matrix.mean(axis=0)
    std = expr_matrix.std(axis=0, ddof=1)
    # 处理标准差为 0 的情况，避免除零
    std = std.where(std != 0, 1.0)
    return (expr_matrix - mean) / std
