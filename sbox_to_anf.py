import random

def truth_table_to_anf_bool(n, truth_table):
    """
    将真值表转换为布尔代数正规型系数（mod 2）。
    
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


def truth_table_to_anf_int(n, truth_table):
    """
    将真值表转换为整数形式ANF系数（整数莫比乌斯变换，不取模2）。

    返回的系数满足：
        f(x) = sum_{S subseteq [n]} a_S * prod_{i in S} x_i
    其中 a_S 为整数。
    """
    m = 1 << n
    coeff = [int(v) for v in truth_table]

    # 整数莫比乌斯变换：若 j 包含第 i 位，则减去去掉该位的系数
    for i in range(n):
        bit = 1 << i
        for j in range(m):
            if j & bit:
                coeff[j] -= coeff[j ^ bit]

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

def current_code_order_masks(n):
    """
    返回当前代码使用的单项式顺序。

    当前顺序为“指数对应整数大小顺序”（mask升序）：
    - 常数项对应mask=0
    - 例如 x2∧x0 对应 mask=4+1=5
    - 例如 x2∧x1∧x0 对应 mask=4+2+1=7
    """
    return list(range(1 << n))


def mask_to_term(mask, n, var_names=None):
    """将mask转换为单项式字符串。mask=0表示常数项1。"""
    if var_names is None:
        var_names = [f"x{i}" for i in range(n)]

    vars_in_term = []
    for j in range(n):
        if mask & (1 << j):
            vars_in_term.append(var_names[j])

    if vars_in_term:
        return " ∧ ".join(vars_in_term)
    return "1"


def anf_expression_bool(n, coeff, var_names=None, masks_order=None):
    """将布尔ANF系数转换为可读的符号表达式（⊕ / ∧）。"""
    if masks_order is None:
        masks_order = current_code_order_masks(n)
    
    terms = []
    for i in masks_order:
        c = coeff[i] & 1
        if c:
            terms.append(mask_to_term(i, n, var_names))
    
    return " ⊕ ".join(terms) if terms else "0"


def anf_expression_int(n, coeff, var_names=None, masks_order=None):
    """将整数ANF系数转换为可读的符号表达式（+ / * / ∧）。"""
    if masks_order is None:
        masks_order = current_code_order_masks(n)

    pieces = []
    for i in masks_order:
        c = coeff[i]
        if c == 0:
            continue

        term = mask_to_term(i, n, var_names)
        if i == 0:
            pieces.append(str(c))
        else:
            if c == 1:
                pieces.append(term)
            elif c == -1:
                pieces.append(f"-{term}")
            else:
                pieces.append(f"{c}*({term})")

    if not pieces:
        return "0"

    expr = " + ".join(pieces)
    return expr.replace("+ -", "- ")
    
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
        
        # 为每个输出位计算两类ANF系数
        # 1) 布尔形式（mod 2）
        # 2) 整数形式（整数莫比乌斯）
        self.anf_coeffs_bool, self.anf_coeffs_int = self._compute_all_anf_coeffs()
        
    def _compute_all_anf_coeffs(self):
        """为S盒的每个输出位计算布尔ANF和整数ANF系数。"""
        anf_coeffs_bool = []
        anf_coeffs_int = []
        
        for output_bit in range(self.m):
            # 提取该输出位的真值表
            truth_table = []
            for input_val in range(self.input_size):
                # 获取S盒输出
                output_val = self.sbox[input_val]
                # 提取第output_bit位
                bit_value = (output_val >> output_bit) & 1
                truth_table.append(bit_value)
            
            # 计算两种ANF系数
            coeff_bool = truth_table_to_anf_bool(self.n, truth_table)
            coeff_int = truth_table_to_anf_int(self.n, truth_table)

            anf_coeffs_bool.append(coeff_bool)
            anf_coeffs_int.append(coeff_int)
        
        return anf_coeffs_bool, anf_coeffs_int
    
    def get_bool_anf_for_bit(self, output_bit):
        """获取指定输出位的布尔ANF系数（mod 2）。"""
        if 0 <= output_bit < self.m:
            return self.anf_coeffs_bool[output_bit]
        else:
            raise ValueError(f"输出位必须在0到{self.m-1}之间")

    def get_int_anf_for_bit(self, output_bit):
        """获取指定输出位的整数ANF系数（整数莫比乌斯）。"""
        if 0 <= output_bit < self.m:
            return self.anf_coeffs_int[output_bit]
        else:
            raise ValueError(f"输出位必须在0到{self.m-1}之间")
    
    def get_bool_anf_expression_for_bit(self, output_bit, var_names=None):
        """获取指定输出位的布尔ANF表达式。"""
        coeff = self.get_bool_anf_for_bit(output_bit)
        return anf_expression_bool(self.n, coeff, var_names)

    def get_int_anf_expression_for_bit(self, output_bit, var_names=None):
        """获取指定输出位的整数ANF表达式。"""
        coeff = self.get_int_anf_for_bit(output_bit)
        return anf_expression_int(self.n, coeff, var_names)

    def get_ordered_coeff_list_for_bit(self, output_bit, mode="bool"):
        """
        获取指定输出位按当前代码顺序（mask整数升序）的全系数列表（含0）。

        mode:
            - "bool": 布尔ANF系数
            - "int": 整数ANF系数
        """
        masks = current_code_order_masks(self.n)
        if mode == "bool":
            coeff = self.get_bool_anf_for_bit(output_bit)
            return [coeff[m] & 1 for m in masks]
        if mode == "int":
            coeff = self.get_int_anf_for_bit(output_bit)
            return [coeff[m] for m in masks]
        raise ValueError("mode 必须是 'bool' 或 'int'")

    def get_ordered_coeff_pairs_for_bit(self, output_bit, mode="bool", var_names=None):
        """返回[(mask, term, coeff)]，用于完整打印每个单项式（含0）。"""
        masks = current_code_order_masks(self.n)
        if mode == "bool":
            coeff = self.get_bool_anf_for_bit(output_bit)
            return [(m, mask_to_term(m, self.n, var_names), coeff[m] & 1) for m in masks]
        if mode == "int":
            coeff = self.get_int_anf_for_bit(output_bit)
            return [(m, mask_to_term(m, self.n, var_names), coeff[m]) for m in masks]
        raise ValueError("mode 必须是 'bool' 或 'int'")
    
    def get_anf_code_for_bit(self, output_bit, output_var=None, var_names=None):
        """获取指定输出位的ANF代码形式"""
        coeff = self.get_bool_anf_for_bit(output_bit)
        
        if output_var is None:
            output_var = f"y{output_bit}"
        
        return anf_to_code(self.n, coeff, output_var, var_names)
    
    def get_all_anf_expressions(self, var_names=None):
        """获取所有输出位的布尔ANF表达式（兼容旧接口）。"""
        expressions = []
        for output_bit in range(self.m):
            expr = self.get_bool_anf_expression_for_bit(output_bit, var_names)
            expressions.append(expr)
        return expressions

    def get_all_bool_anf_expressions(self, var_names=None):
        expressions = []
        for output_bit in range(self.m):
            expressions.append(self.get_bool_anf_expression_for_bit(output_bit, var_names))
        return expressions

    def get_all_int_anf_expressions(self, var_names=None):
        expressions = []
        for output_bit in range(self.m):
            expressions.append(self.get_int_anf_expression_for_bit(output_bit, var_names))
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

    def print_analysis(self, var_names=None, show_code=True, show_full_coeffs=True):
        """打印完整的S盒分析结果"""
        print(f"=== S盒 (n={self.n}比特输入, m={self.m}比特输出) ===")
        #print(f"S盒值 (前16个): {self.sbox[:16]}{'...' if len(self.sbox) > 16 else ''}")
        print()
        
        print("输出位的布尔ANF表达式（mod2）:")
        expressions_bool = self.get_all_bool_anf_expressions(var_names)
        for i, expr in enumerate(expressions_bool):
            print(f"  输出位 y{i}: {expr}")

        print()
        print("输出位的整数ANF表达式（整数莫比乌斯）:")
        expressions_int = self.get_all_int_anf_expressions(var_names)
        for i, expr in enumerate(expressions_int):
            print(f"  输出位 y{i}: {expr}")
        
        print()

        if show_full_coeffs:
            print("按当前代码顺序（指数对应整数大小，即mask升序）输出全系数（含0）:")
            for output_bit in range(self.m):
                print(f"  // 输出位 y{output_bit} 布尔ANF系数（0..255）")
                print(f"  {self.get_ordered_coeff_list_for_bit(output_bit, mode='bool')}")
                print(f"  // 输出位 y{output_bit} 整数ANF系数（0..255）")
                print(f"  {self.get_ordered_coeff_list_for_bit(output_bit, mode='int')}")
                print()

            print("按当前代码顺序逐项打印（mask, term, coeff），含0:")
            for output_bit in range(self.m):
                print(f"  // 输出位 y{output_bit} 整数ANF逐项系数")
                pairs = self.get_ordered_coeff_pairs_for_bit(output_bit, mode="int", var_names=var_names)
                for mask, term, coeff in pairs:
                    print(f"  {mask:3d}: {coeff:4d}  // {term}")
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
    SM4-sbox = [0xd6, 0x90, 0xe9, 0xfe, 0xcc, 0xe1, 0x3d, 0xb7, 0x16, 0xb6, 0x14, 0xc2, 0x28, 0xfb, 0x2c, 0x05,
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
    AES-sbox = [0x63, 0x7C, 0x77, 0x7B, 0xF2, 0x6B, 0x6F, 0xC5,
            0x30, 0x01, 0x67, 0x2B, 0xFE, 0xD7, 0xAB, 0x76,
            0xCA, 0x82, 0xC9, 0x7D, 0xFA, 0x59, 0x47, 0xF0,
            0xAD, 0xD4, 0xA2, 0xAF, 0x9C, 0xA4, 0x72, 0xC0,
            0xB7, 0xFD, 0x93, 0x26, 0x36, 0x3F, 0xF7, 0xCC,
            0x34, 0xA5, 0xE5, 0xF1, 0x71, 0xD8, 0x31, 0x15,
            0x04, 0xC7, 0x23, 0xC3, 0x18, 0x96, 0x05, 0x9A,
            0x07, 0x12, 0x80, 0xE2, 0xEB, 0x27, 0xB2, 0x75,
            0x09, 0x83, 0x2C, 0x1A, 0x1B, 0x6E, 0x5A, 0xA0,
            0x52, 0x3B, 0xD6, 0xB3, 0x29, 0xE3, 0x2F, 0x84,
            0x53, 0xD1, 0x00, 0xED, 0x20, 0xFC, 0xB1, 0x5B,
            0x6A, 0xCB, 0xBE, 0x39, 0x4A, 0x4C, 0x58, 0xCF,
            0xD0, 0xEF, 0xAA, 0xFB, 0x43, 0x4D, 0x33, 0x85,
            0x45, 0xF9, 0x02, 0x7F, 0x50, 0x3C, 0x9F, 0xA8,
            0x51, 0xA3, 0x40, 0x8F, 0x92, 0x9D, 0x38, 0xF5,
            0xBC, 0xB6, 0xDA, 0x21, 0x10, 0xFF, 0xF3, 0xD2,
            0xCD, 0x0C, 0x13, 0xEC, 0x5F, 0x97, 0x44, 0x17,
            0xC4, 0xA7, 0x7E, 0x3D, 0x64, 0x5D, 0x19, 0x73,
            0x60, 0x81, 0x4F, 0xDC, 0x22, 0x2A, 0x90, 0x88,
            0x46, 0xEE, 0xB8, 0x14, 0xDE, 0x5E, 0x0B, 0xDB,
            0xE0, 0x32, 0x3A, 0x0A, 0x49, 0x06, 0x24, 0x5C,
            0xC2, 0xD3, 0xAC, 0x62, 0x91, 0x95, 0xE4, 0x79,
            0xE7, 0xC8, 0x37, 0x6D, 0x8D, 0xD5, 0x4E, 0xA9,
            0x6C, 0x56, 0xF4, 0xEA, 0x65, 0x7A, 0xAE, 0x08,
            0xBA, 0x78, 0x25, 0x2E, 0x1C, 0xA6, 0xB4, 0xC6,
            0xE8, 0xDD, 0x74, 0x1F, 0x4B, 0xBD, 0x8B, 0x8A,
            0x70, 0x3E, 0xB5, 0x66, 0x48, 0x03, 0xF6, 0x0E,
            0x61, 0x35, 0x57, 0xB9, 0x86, 0xC1, 0x1D, 0x9E,
            0xE1, 0xF8, 0x98, 0x11, 0x69, 0xD9, 0x8E, 0x94,
            0x9B, 0x1E, 0x87, 0xE9, 0xCE, 0x55, 0x28, 0xDF,
            0x8C, 0xA1, 0x89, 0x0D, 0xBF, 0xE6, 0x42, 0x68,
            0x41, 0x99, 0x2D, 0x0F, 0xB0, 0x54, 0xBB, 0x16];
    
    analyzer = SBoxAnalyzer(n, m, AES-sbox)
    analyzer.print_analysis()
    
    print("\n" + "=" * 60)
