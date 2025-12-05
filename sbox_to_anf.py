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
    
def anf_to_code(n, coeff, output_var="y", var_names=None):
    """将ANF系数转换为可执行的代码形式"""
    if var_names is None:
        var_names = [f"x{i}" for i in range(n)]
    
    lines = []
    
    # 初始化输出变量
    lines.append(f"{output_var} = 0;")
    
    # 处理常数项1
    if coeff[0]:
        lines.append(f"{output_var} ^= 1;")
    
    # 处理所有非常数项
    for i in range(1, 1 << n):
        if coeff[i]:
            # 获取变量索引
            vars_in_term = []
            for j in range(n):
                if i & (1 << j):
                    vars_in_term.append(var_names[j])
            
            if len(vars_in_term) == 1:
                # 单个变量
                lines.append(f"{output_var} ^= {vars_in_term[0]};")
            else:
                # 多个变量的与运算
                and_expr = " & ".join(vars_in_term)
                lines.append(f"{output_var} ^= ({and_expr});")
    
    # 如果没有任何项，输出为0
    if len(lines) == 1 and not coeff[0]:
        lines[0] = f"{output_var} = 0;"
    
    return "\n".join(lines)

def anf_expression_to_code(expr, output_var="y", var_names=None):
    """将ANF表达式字符串转换为可执行的代码形式"""
    # 从表达式解析项
    terms = expr.split(" ⊕ ")
    lines = []
    
    lines.append(f"{output_var} = 0;")
    
    for term in terms:
        if term == "1":
            lines.append(f"{output_var} ^= 1;")
        elif term == "0":
            # 常数0，什么都不做
            continue
        elif " ∧ " in term:
            # 多个变量的与
            vars_in_term = term.split(" ∧ ")
            and_expr = " & ".join(vars_in_term)
            lines.append(f"{output_var} ^= ({and_expr});")
        else:
            # 单个变量
            lines.append(f"{output_var} ^= {term};")
    
    return "\n".join(lines)
    
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
    
    def get_anf_code_for_bit(self, output_bit, output_var=None, var_names=None):
        """获取指定输出位的ANF代码形式"""
        coeff = self.get_anf_for_bit(output_bit)
        
        if output_var is None:
            output_var = f"y{output_bit}"
        
        return anf_to_code(self.n, coeff, output_var, var_names)
    
    def get_all_anf_expressions(self, var_names=None):
        """获取所有输出位的ANF表达式"""
        expressions = []
        for output_bit in range(self.m):
            expr = self.get_anf_expression_for_bit(output_bit, var_names)
            expressions.append(expr)
        return expressions
    
    def get_all_anf_codes(self, output_vars=None, var_names=None):
        """获取所有输出位的ANF代码形式"""
        codes = []
        for output_bit in range(self.m):
            if output_vars is not None and output_bit < len(output_vars):
                output_var = output_vars[output_bit]
            else:
                output_var = f"y{output_bit}"
            
            code = self.get_anf_code_for_bit(output_bit, output_var, var_names)
            codes.append(code)
        return codes

    def print_analysis(self, var_names=None, show_code=True):
        """打印完整的S盒分析结果"""
        print(f"=== S盒 (n={self.n}比特输入, m={self.m}比特输出) ===")
        #print(f"S盒值 (前16个): {self.sbox[:16]}{'...' if len(self.sbox) > 16 else ''}")
        print()
        
        print("输出位的ANF表达式:")
        expressions = self.get_all_anf_expressions(var_names)
        for i, expr in enumerate(expressions):
            print(f"  输出位 y{i}: {expr}")
        
        print()
        
        if show_code:
            print("ANF代码形式:")
            codes = self.get_all_anf_codes(var_names=var_names)
            for i, code in enumerate(codes):
                print(f"  // 输出位 y{i}")
                for line in code.split('\n'):
                    print(f"  {line}")
                print()

if __name__ == "__main__":
    
    print("=" * 60)
    
    n, m = 8, 8
    sbox = [0xd6, 0x90, 0xe9, 0xfe, 0xcc, 0xe1, 0x3d, 0xb7, 0x16, 0xb6, 0x14, 0xc2, 0x28, 0xfb, 0x2c, 0x05,
    0x2b, 0x67, 0x9a, 0x76, 0x2a, 0xbe, 0x04, 0xc3, 0xaa, 0x44, 0x13, 0x26, 0x49, 0x86, 0x06, 0x99,
    0x9c, 0x42, 0x50, 0xf4, 0x91, 0xef, 0x98, 0x7a, 0x33, 0x54, 0x0b, 0x43, 0xed, 0xcf, 0xac, 0x62,
    0xe4, 0xb3, 0x1c, 0xa9, 0xc9, 0x08, 0xe8, 0x95, 0x80, 0xdf, 0x94, 0xfa, 0x75, 0x8f, 0x3f, 0xa6,
    0x47, 0x07, 0xa7, 0xfc, 0xf3, 0x73, 0x17, 0xba, 0x83, 0x59, 0x3c, 0x19, 0xe6, 0x85, 0x4f, 0xa8,
    0x68, 0x6b, 0x81, 0xb2, 0x71, 0x64, 0xda, 0x8b, 0xf8, 0xeb, 0x0f, 0x4b, 0x70, 0x56, 0x9d, 0x35,
    0x1e, 0x24, 0x0e, 0x5e, 0x63, 0x58, 0xd1, 0xa2, 0x25, 0x22, 0x7c, 0x3b, 0x01, 0x21, 0x78, 0x87,
    0xd4, 0x00, 0x46, 0x57, 0x9f, 0xd3, 0x27, 0x52, 0x4c, 0x36, 0x02, 0xe7, 0xa0, 0xc4, 0xc8, 0x9e,
    0xea, 0xbf, 0x8a, 0xd2, 0x40, 0xc7, 0x38, 0xb5, 0xa3, 0xf7, 0xf2, 0xce, 0xf9, 0x61, 0x15, 0xa1,
    0xe0, 0xae, 0x5d, 0xa4, 0x9b, 0x34, 0x1a, 0x55, 0xad, 0x93, 0x32, 0x30, 0xf5, 0x8c, 0xb1, 0xe3,
    0x1d, 0xf6, 0xe2, 0x2e, 0x82, 0x66, 0xca, 0x60, 0xc0, 0x29, 0x23, 0xab, 0x0d, 0x53, 0x4e, 0x6f,
    0xd5, 0xdb, 0x37, 0x45, 0xde, 0xfd, 0x8e, 0x2f, 0x03, 0xff, 0x6a, 0x72, 0x6d, 0x6c, 0x5b, 0x51,
    0x8d, 0x1b, 0xaf, 0x92, 0xbb, 0xdd, 0xbc, 0x7f, 0x11, 0xd9, 0x5c, 0x41, 0x1f, 0x10, 0x5a, 0xd8,
    0x0a, 0xc1, 0x31, 0x88, 0xa5, 0xcd, 0x7b, 0xbd, 0x2d, 0x74, 0xd0, 0x12, 0xb8, 0xe5, 0xb4, 0xb0,
    0x89, 0x69, 0x97, 0x4a, 0x0c, 0x96, 0x77, 0x7e, 0x65, 0xb9, 0xf1, 0x09, 0xc5, 0x6e, 0xc6, 0x84,
    0x18, 0xf0, 0x7d, 0xec, 0x3a, 0xdc, 0x4d, 0x20, 0x79, 0xee, 0x5f, 0x3e, 0xd7, 0xcb, 0x39, 0x48];
    
    analyzer = SBoxAnalyzer(n, m, sbox)
    analyzer.print_analysis()
    
    print("\n" + "=" * 60)
