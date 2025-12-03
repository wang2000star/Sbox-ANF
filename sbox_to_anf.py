import random
import numpy as np

def truth_table_to_anf(n, truth_table):
    """
    将真值表转换为代数正规型系数。
    
    参数:
        n: 变量个数
        truth_table: 列表，长度为 2^n，按自然二进制顺序排列的真值表输出
        
    返回:
        coeff: ANF系数列表，长度为 2^n
    """
    m = 1 << n
    coeff = truth_table.copy()
    
    # 快速莫比乌斯变换
    for i in range(n):
        mask = 1 << i
        for j in range(m):
            if j & mask:  # 只处理第i位为1的位置
                coeff[j] ^= coeff[j ^ mask]
    
    return coeff

def anf_to_truth_table(n, coeff):
    """从ANF系数重建真值表"""
    m = 1 << n
    table = [0] * m
    for i in range(m):
        # 计算所有满足 j ⊆ i 的系数异或
        j = i
        while True:
            table[i] ^= coeff[j]
            if j == 0:
                break
            j = (j - 1) & i
    return table

def anf_expression(n, coeff, var_names=None):
    """将ANF系数转换为可读的表达式"""
    if var_names is None:
        var_names = [f"x{i}" for i in range(n)]
    
    terms = []
    for i in range(1 << n):
        if coeff[i]:
            # 获取变量名
            vars_in_term = []
            for j in range(n):
                if i & (1 << j):
                    vars_in_term.append(var_names[j])
            if vars_in_term:
                term = " ∧ ".join(vars_in_term)
                terms.append(term)
            else:
                terms.append("1")  # 常数项1
    
    return " ⊕ ".join(terms) if terms else "0"

class SBoxAnalyzer:
    """S盒分析器，分析n比特输入m比特输出的S盒"""
    
    def __init__(self, n, m, sbox=None):
        """
        初始化S盒分析器
        
        参数:
            n: 输入比特数
            m: 输出比特数
            sbox: S盒映射，长度为2^n的列表，每个元素是0到2^m-1之间的整数
        """
        self.n = n
        self.m = m
        self.input_size = 1 << n
        self.output_size = 1 << m
        
        if sbox is None:
            # 如果没有提供S盒，生成一个随机映射（不是置换）
            self.sbox = [random.randint(0, self.output_size - 1) for _ in range(self.input_size)]
        else:
            self.sbox = sbox.copy()
        
        # 为每个输出位计算ANF系数
        self.anf_coeffs = self._compute_all_anf_coeffs()
        
    def _compute_all_anf_coeffs(self):
        """为S盒的每个输出位计算ANF系数"""
        anf_coeffs = []
        
        for output_bit in range(self.m):
            # 提取该输出位的真值表
            truth_table = []
            for input_val in range(self.input_size):
                # 获取S盒输出
                output_val = self.sbox[input_val]
                # 提取第output_bit位
                bit_value = (output_val >> output_bit) & 1
                truth_table.append(bit_value)
            
            # 计算ANF系数
            coeff = truth_table_to_anf(self.n, truth_table)
            anf_coeffs.append(coeff)
        
        return anf_coeffs
    
    def get_anf_for_bit(self, output_bit):
        """获取指定输出位的ANF系数"""
        if 0 <= output_bit < self.m:
            return self.anf_coeffs[output_bit]
        else:
            raise ValueError(f"输出位必须在0到{self.m-1}之间")
    
    def get_anf_expression_for_bit(self, output_bit, var_names=None):
        """获取指定输出位的ANF表达式"""
        coeff = self.get_anf_for_bit(output_bit)
        return anf_expression(self.n, coeff, var_names)
    
    def get_all_anf_expressions(self, var_names=None):
        """获取所有输出位的ANF表达式"""
        expressions = []
        for output_bit in range(self.m):
            expr = self.get_anf_expression_for_bit(output_bit, var_names)
            expressions.append(expr)
        return expressions
    
    def algebraic_degree_for_bit(self, output_bit):
        """计算指定输出位的代数次数"""
        coeff = self.get_anf_for_bit(output_bit)
        
        max_degree = 0
        for i in range(self.input_size):
            if coeff[i]:
                # 计算单项式的次数（i的二进制表示中1的个数）
                degree = bin(i).count('1')
                if degree > max_degree:
                    max_degree = degree
        
        return max_degree
    
    def algebraic_degrees(self):
        """计算所有输出位的代数次数"""
        degrees = []
        for output_bit in range(self.m):
            degree = self.algebraic_degree_for_bit(output_bit)
            degrees.append(degree)
        return degrees
    
    def nonlinearity_for_bit(self, output_bit):
        """计算指定输出位的非线性度（到所有仿射函数的汉明距离的最小值）"""
        coeff = self.get_anf_for_bit(output_bit)
        truth_table = anf_to_truth_table(self.n, coeff)
        
        min_distance = self.input_size  # 初始化为最大值
        
        # 检查所有仿射函数：c0 ⊕ c1*x1 ⊕ ... ⊕ cn*xn
        # 共有 2^(n+1) 个仿射函数
        for a in range(1 << (self.n + 1)):
            # a的二进制表示：最低位是c0，接下来的n位是c1到cn
            affine_table = []
            for x in range(self.input_size):
                # 计算仿射函数值
                value = (a & 1)  # c0
                for i in range(self.n):
                    if (a >> (i + 1)) & 1:  # 如果ci=1
                        # 提取xi的值
                        xi = (x >> i) & 1
                        value ^= xi
                affine_table.append(value)
            
            # 计算汉明距离
            distance = sum(1 for i in range(self.input_size) if truth_table[i] != affine_table[i])
            
            if distance < min_distance:
                min_distance = distance
        
        return min_distance
    
    def nonlinearities(self):
        """计算所有输出位的非线性度"""
        nonlinearities = []
        for output_bit in range(self.m):
            nl = self.nonlinearity_for_bit(output_bit)
            nonlinearities.append(nl)
        return nonlinearities
    
    def print_analysis(self, var_names=None):
        """打印完整的S盒分析结果"""
        print(f"=== S盒分析 (n={self.n}比特输入, m={self.m}比特输出) ===")
        print(f"S盒值 (前16个): {self.sbox[:16]}{'...' if len(self.sbox) > 16 else ''}")
        print()
        
        print("输出位的ANF表达式:")
        expressions = self.get_all_anf_expressions(var_names)
        for i, expr in enumerate(expressions):
            print(f"  输出位 y{i}: {expr}")
        
        print()
        
        print("代数次数:")
        degrees = self.algebraic_degrees()
        for i, degree in enumerate(degrees):
            print(f"  输出位 y{i}: {degree}")
        print(f"  平均代数次数: {sum(degrees)/len(degrees):.2f}")
        
        print()
        
        print("非线性度:")
        nls = self.nonlinearities()
        for i, nl in enumerate(nls):
            print(f"  输出位 y{i}: {nl}")
        print(f"  平均非线性度: {sum(nls)/len(nls):.2f}")
        
        print()
        
        # 打印真值表示例
        print("真值表示例 (前8个输入):")
        print("输入\t二进制\t\tS盒输出\t\t二进制")
        for i in range(min(8, self.input_size)):
            sbox_output = self.sbox[i]
            input_bin = format(i, f'0{self.n}b')
            output_bin = format(sbox_output, f'0{self.m}b')
            print(f"{i}\t{input_bin}\t\t{sbox_output}\t\t{output_bin}")

def test_specific_sboxes():
    """测试特定的S盒"""
    
    print("=" * 60)
    print("测试1: 4x4 S盒 (n=4, m=4)")
    print("=" * 60)
    
    # 示例1: 简单的4x4 S盒
    n, m = 4, 4
    
    # 创建一个简单的S盒（这里使用一个随机置换）
    sbox = list(range(16))
    random.shuffle(sbox)
    
    analyzer = SBoxAnalyzer(n, m, sbox)
    analyzer.print_analysis()
    
    print("\n" + "=" * 60)
    print("测试2: 8x8 S盒 (AES的S盒类似)")
    print("=" * 60)
    
    # 示例2: 8x8 S盒（仿照AES的S盒，这里用随机值代替）
    n, m = 8, 8
    
    # 生成一个伪随机的8x8 S盒
    sbox8 = list(range(256))
    random.shuffle(sbox8)
    
    analyzer8 = SBoxAnalyzer(n, m, sbox8)
    
    # 对于8x8 S盒，只显示前几个输出位的ANF表达式
    print(f"S盒大小: {len(sbox8)}")
    print("\n前3个输出位的ANF表达式:")
    for i in range(3):
        expr = analyzer8.get_anf_expression_for_bit(i)
        print(f"  输出位 y{i}: {expr}")
    
    print("\n所有输出位的代数次数:")
    degrees = analyzer8.algebraic_degrees()
    for i, degree in enumerate(degrees):
        print(f"  y{i}: {degree}", end="  ")
        if (i + 1) % 8 == 0:
            print()
    
    print(f"\n平均代数次数: {sum(degrees)/len(degrees):.2f}")
    
    print("\n所有输出位的非线性度:")
    nls = analyzer8.nonlinearities()
    for i, nl in enumerate(nls):
        print(f"  y{i}: {nl}", end="  ")
        if (i + 1) % 8 == 0:
            print()
    
    print(f"\n平均非线性度: {sum(nls)/len(nls):.2f}")
    
    print("\n" + "=" * 60)
    print("测试3: 非方形的S盒 (6x4 S盒)")
    print("=" * 60)
    
    # 示例3: 非方形S盒
    n, m = 6, 4
    sbox_6x4 = [random.randint(0, 15) for _ in range(64)]
    
    analyzer_6x4 = SBoxAnalyzer(n, m, sbox_6x4)
    
    print("所有输出位的ANF表达式:")
    expressions = analyzer_6x4.get_all_anf_expressions()
    for i, expr in enumerate(expressions):
        print(f"  输出位 y{i}: {expr}")
    
    print("\n代数次数:", analyzer_6x4.algebraic_degrees())
    print("非线性度:", analyzer_6x4.nonlinearities())

def test_affine_function():
    """测试仿射函数作为S盒的情况"""
    print("\n" + "=" * 60)
    print("测试: 仿射函数S盒")
    print("=" * 60)
    
    # 创建一个仿射函数: y = x0 ⊕ x1 ⊕ x2 ⊕ 1 (3比特输入，1比特输出)
    n, m = 3, 1
    sbox_affine = []
    for i in range(8):
        # 提取输入位
        x0 = (i >> 0) & 1
        x1 = (i >> 1) & 1
        x2 = (i >> 2) & 1
        # 仿射函数: y = x0 ⊕ x1 ⊕ x2 ⊕ 1
        y = x0 ^ x1 ^ x2 ^ 1
        sbox_affine.append(y)
    
    analyzer = SBoxAnalyzer(n, m, sbox_affine)
    
    print("S盒值:", sbox_affine)
    print("ANF表达式:", analyzer.get_anf_expression_for_bit(0))
    print("代数次数:", analyzer.algebraic_degree_for_bit(0))
    print("非线性度:", analyzer.nonlinearity_for_bit(0))
    
    # 验证ANF表达式是否正确
    expected_expr = "1 ⊕ x0 ⊕ x1 ⊕ x2"
    actual_expr = analyzer.get_anf_expression_for_bit(0)
    print(f"表达式正确: {actual_expr == expected_expr}")

def test_complete_analysis():
    """完整的S盒分析测试"""
    
    # 创建一个自定义的4x4 S盒
    n, m = 4, 4
    
    # 定义一个S盒（这里使用PRESENT密码算法的S盒）
    present_sbox = [
        0xC, 0x5, 0x6, 0xB, 0x9, 0x0, 0xA, 0xD,
        0x3, 0xE, 0xF, 0x8, 0x4, 0x7, 0x1, 0x2
    ]
    
    print("=" * 60)
    print("PRESENT密码算法S盒分析")
    print("=" * 60)
    
    analyzer = SBoxAnalyzer(n, m, present_sbox)
    
    # 获取所有分析结果
    expressions = analyzer.get_all_anf_expressions()
    degrees = analyzer.algebraic_degrees()
    nls = analyzer.nonlinearities()
    
    print("S盒值:", present_sbox)
    print("\n二进制表示:")
    for i, val in enumerate(present_sbox):
        print(f"  S({i:2d}) = {val:2d} = {format(val, '04b')}")
    
    print("\nANF表达式:")
    for i, expr in enumerate(expressions):
        print(f"  y{i}: {expr}")
    
    print(f"\n代数次数: {degrees}")
    print(f"非线性度: {nls}")
    
    # 验证函数正确性
    print("\n验证: 从ANF重建S盒输出位")
    for output_bit in range(m):
        coeff = analyzer.get_anf_for_bit(output_bit)
        reconstructed = anf_to_truth_table(n, coeff)
        
        # 提取实际的S盒输出位
        actual = []
        for input_val in range(1 << n):
            output_val = present_sbox[input_val]
            bit_value = (output_val >> output_bit) & 1
            actual.append(bit_value)
        
        if reconstructed == actual:
            print(f"  输出位 y{output_bit}: 正确")
        else:
            print(f"  输出位 y{output_bit}: 错误")

if __name__ == "__main__":
    print("S盒布尔函数分析工具")
    print("=" * 60)
    
    # 运行各种测试
    test_specific_sboxes()
    test_affine_function()
    test_complete_analysis()
    
    # 性能测试：分析一个较大的S盒
    print("\n" + "=" * 60)
    print("性能测试: 随机10x8 S盒")
    print("=" * 60)
    
    n, m = 10, 8
    large_sbox = [random.randint(0, (1 << m) - 1) for _ in range(1 << n)]
    
    import time
    start_time = time.time()
    
    analyzer = SBoxAnalyzer(n, m, large_sbox)
    
    # 只计算前几个输出位的ANF表达式
    print(f"S盒大小: {len(large_sbox)}")
    print("前3个输出位的代数次数:")
    for i in range(3):
        degree = analyzer.algebraic_degree_for_bit(i)
        print(f"  输出位 y{i}: {degree}")
    
    end_time = time.time()
    print(f"\n分析时间: {end_time - start_time:.4f}秒")
