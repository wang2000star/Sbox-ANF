import random

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

def count_anf_terms(coeff):
    """统计ANF中的项数（非零系数个数）"""
    return sum(1 for c in coeff if c)

def test_random_truth_table(n, num_tests=10):
    """测试随机生成的真值表"""
    print(f"===== 测试 {n} 变量布尔函数 =====")
    print(f"每个测试生成 {1<<n} 个真值表输出")
    print("=" * 50)
    
    for test_num in range(1, num_tests + 1):
        # 生成随机真值表
        m = 1 << n
        truth_table = [random.randint(0, 1) for _ in range(m)]
        
        # 计算ANF系数
        coeff = truth_table_to_anf(n, truth_table)
        
        # 重建真值表
        reconstructed = anf_to_truth_table(n, coeff)
        
        # 验证一致性
        is_correct = reconstructed == truth_table
        
        # 统计信息
        num_terms = count_anf_terms(coeff)
        
        print(f"测试 #{test_num}:")
        print(f"  真值表随机生成 - 前10个值: {truth_table[:10]}...")
        print(f"  ANF系数项数: {num_terms}/{m} (非零系数/总系数)")
        print(f"  转换正确: {is_correct}")
        
        if not is_correct:
            print("  ERROR: 转换错误!")
            return False
        
        # 显示部分ANF项
        if num_terms <= 10:  # 如果项数不多，显示所有项
            terms = []
            for i in range(m):
                if coeff[i]:
                    var_terms = [f"x{j}" for j in range(n) if i & (1 << j)]
                    term = " ∧ ".join(var_terms) if var_terms else "1"
                    terms.append(term)
            print(f"  ANF表达式: {' ⊕ '.join(terms)}")
        else:
            # 只显示前5项
            terms_shown = 0
            term_list = []
            for i in range(m):
                if coeff[i]:
                    if terms_shown < 5:
                        var_terms = [f"x{j}" for j in range(n) if i & (1 << j)]
                        term = " ∧ ".join(var_terms) if var_terms else "1"
                        term_list.append(term)
                        terms_shown += 1
                    else:
                        break
            print(f"  ANF表达式 (前5项): {' ⊕ '.join(term_list)} ⊕ ...")
        
        print("-" * 50)
    
    print(f"所有 {num_tests} 个测试都通过!")
    return True

def test_specific_functions():
    """测试特定的布尔函数"""
    print("===== 测试特定函数 =====")
    
    # 1. 常数0函数
    n = 4
    truth_table = [0] * (1 << n)
    coeff = truth_table_to_anf(n, truth_table)
    print(f"常数0函数: ANF项数={count_anf_terms(coeff)}")
    
    # 2. 常数1函数
    truth_table = [1] * (1 << n)
    coeff = truth_table_to_anf(n, truth_table)
    print(f"常数1函数: ANF项数={count_anf_terms(coeff)}")
    
    # 3. 奇偶校验函数
    n = 4
    truth_table = []
    for i in range(1 << n):
        # 计算i的二进制表示中1的个数
        parity = bin(i).count('1') % 2
        truth_table.append(parity)
    coeff = truth_table_to_anf(n, truth_table)
    print(f"{n}变量奇偶校验函数: ANF项数={count_anf_terms(coeff)}")
    
    # 4. AND函数 (所有变量的与)
    n = 4
    truth_table = []
    for i in range(1 << n):
        # 只有当所有位都是1时才输出1
        truth_table.append(1 if i == (1 << n) - 1 else 0)
    coeff = truth_table_to_anf(n, truth_table)
    print(f"{n}变量AND函数: ANF项数={count_anf_terms(coeff)}")
    
    # 5. OR函数 (所有变量的或)
    n = 4
    truth_table = []
    for i in range(1 << n):
        # 当至少有一位是1时输出1
        truth_table.append(1 if i > 0 else 0)
    coeff = truth_table_to_anf(n, truth_table)
    print(f"{n}变量OR函数: ANF项数={count_anf_terms(coeff)}")

def performance_test(n=8):
    """性能测试"""
    print("===== 性能测试 =====")
    m = 1 << n
    print(f"测试 {n} 变量函数 ({m} 个真值表输出)")
    
    # 生成随机真值表
    truth_table = [random.randint(0, 1) for _ in range(m)]
    
    # 测试转换时间
    import time
    start_time = time.time()
    coeff = truth_table_to_anf(n, truth_table)
    end_time = time.time()
    
    print(f"转换时间: {end_time - start_time:.6f} 秒")
    print(f"ANF系数非零项数: {count_anf_terms(coeff)}/{m}")
    
    # 验证转换正确性
    reconstructed = anf_to_truth_table(n, coeff)
    print(f"转换正确: {reconstructed == truth_table}")

# 运行测试
if __name__ == "__main__":
    print("布尔真值表到ANF表达式转换测试")
    print("=" * 60)
    
    # 测试小规模函数
    test_specific_functions()
    print()
    
    # 测试3变量随机真值表
    test_random_truth_table(3, num_tests=3)
    print()
    
    # 测试4变量随机真值表
    test_random_truth_table(4, num_tests=2)
    print()
    
    # 测试8变量随机真值表 (256个输出)
    test_random_truth_table(8, num_tests=2)
    print()
    
    # 性能测试
    performance_test(8)
    print()
    
    # 测试大规模转换
    print("===== 大规模测试 =====")
    for n in [10, 12]:
        m = 1 << n
        print(f"\n{n} 变量函数 ({m} 个输出):")
        truth_table = [random.randint(0, 1) for _ in range(m)]
        coeff = truth_table_to_anf(n, truth_table)
        reconstructed = anf_to_truth_table(n, coeff)
        print(f"  正确性: {reconstructed == truth_table}")
        print(f"  ANF项数: {count_anf_terms(coeff)}/{m}")
