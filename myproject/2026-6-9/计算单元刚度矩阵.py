import numpy as np

def truss3d_element_stiffness(x1, x2, E, A):
    """
    计算三维杆单元的单元刚度矩阵
    x1, x2: 变形前节点坐标 [x, y, z]
    E: 弹性模量
    A: 截面积
    返回: Ke, L, direction_cosines
    """
    delta = np.array(x2) - np.array(x1)
    L = np.linalg.norm(delta)
    
    if L < 1e-12:
        raise ValueError("错误：两个节点重合，单元退化！")
    
    # 方向余弦
    cx, cy, cz = delta / L
    
    # 刚度矩阵系数
    k = E * A / L
    
    # 构建 6x6 刚度矩阵
    Ke = k * np.array([
        [ cx*cx,  cx*cy,  cx*cz, -cx*cx, -cx*cy, -cx*cz],
        [ cx*cy,  cy*cy,  cy*cz, -cx*cy, -cy*cy, -cy*cz],
        [ cx*cz,  cy*cz,  cz*cz, -cx*cz, -cy*cz, -cz*cz],
        [-cx*cx, -cx*cy, -cx*cz,  cx*cx,  cx*cy,  cx*cz],
        [-cx*cy, -cy*cy, -cy*cz,  cx*cy,  cy*cy,  cy*cz],
        [-cx*cz, -cy*cz, -cz*cz,  cx*cz,  cy*cz,  cz*cz]
    ])
    
    return Ke, L, (cx, cy, cz)


def compute_displacement(x1, x2, x1_prime, x2_prime):
    """
    由变形前后坐标计算节点位移
    x1, x2: 变形前坐标
    x1_prime, x2_prime: 变形后坐标
    返回: de = [u1, v1, w1, u2, v2, w2]
    """
    de1 = np.array(x1_prime) - np.array(x1)
    de2 = np.array(x2_prime) - np.array(x2)
    de = np.concatenate([de1, de2])
    return de


def compute_stress_force_from_coords(x1, x2, x1_prime, x2_prime, E, A):
    """
    基于变形前后坐标直接计算应变、应力、轴力
    (避免先计算位移再应变，更直接)
    """
    # 变形前长度和方向
    delta_orig = np.array(x2) - np.array(x1)
    L_orig = np.linalg.norm(delta_orig)
    
    if L_orig < 1e-12:
        raise ValueError("错误：变形前两个节点重合！")
    
    cx, cy, cz = delta_orig / L_orig
    
    # 变形后长度
    delta_def = np.array(x2_prime) - np.array(x1_prime)
    L_def = np.linalg.norm(delta_def)
    
    # 工程应变 (Green-Lagrange 简化为小变形)
    epsilon = (L_def - L_orig) / L_orig
    
    sigma = E * epsilon
    N = sigma * A
    
    return epsilon, sigma, N, L_orig, L_def, (cx, cy, cz)


def check_stiffness_matrix(Ke):
    """
    检查刚度矩阵的性质
    """
    checks = {}
    
    # 1. 对称性
    checks['对称性'] = np.allclose(Ke, Ke.T, atol=1e-10)
    
    # 2. 对角元素非负性
    diag = np.diag(Ke)
    checks['对角元素非负'] = np.all(diag >= -1e-10)
    checks['对角元素值'] = diag
    
    # 3. 特征值
    eigvals = np.linalg.eigvalsh(Ke)
    checks['特征值'] = eigvals
    checks['零特征值个数'] = np.sum(np.abs(eigvals) < 1e-10)
    checks['最小特征值'] = np.min(eigvals)
    
    # 4. 行列式
    det = np.linalg.det(Ke)
    checks['行列式'] = det
    checks['是否奇异'] = np.abs(det) < 1e-10
    
    # 5. 秩
    rank = np.linalg.matrix_rank(Ke, tol=1e-10)
    checks['秩'] = rank
    checks['秩亏损'] = 6 - rank
    
    return checks


# ===================== 主程序 =====================
if __name__ == "__main__":
    print("=" * 80)
    print("三维杆单元计算程序 (基于变形前后坐标)")
    print("=" * 80)
    
    # ========== 输入变形前坐标 ==========
    print("\n【1. 输入变形前节点坐标】")
    coords_orig = list(map(float, input(
        "变形前坐标 [x1, y1, z1, x2, y2, z2] 用空格分隔: ").split()))
    
    if len(coords_orig) != 6:
        raise ValueError("需要输入6个坐标值 (x1,y1,z1,x2,y2,z2)！")
    
    x1 = coords_orig[0:3]
    x2 = coords_orig[3:6]
    
    # ========== 输入变形后坐标 ==========
    print("\n【2. 输入变形后节点坐标】")
    coords_def = list(map(float, input(
        "变形后坐标 [x1', y1', z1', x2', y2', z2'] 用空格分隔: ").split()))
    
    if len(coords_def) != 6:
        raise ValueError("需要输入6个坐标值 (x1',y1',z1',x2',y2',z2')！")
    
    x1_prime = coords_def[0:3]
    x2_prime = coords_def[3:6]
    
    # ========== 输入材料参数 ==========
    print("\n【3. 输入材料与截面参数】")
    E = float(input("弹性模量 E (Pa): "))
    A = float(input("截面积 A (m^2): "))
    
    # ========== 计算节点位移 ==========
    de = compute_displacement(x1, x2, x1_prime, x2_prime)
    
    print("\n" + "=" * 80)
    print("计算结果")
    print("=" * 80)
    
    # ========== 输出位移 ==========
    print(f"\n【节点位移】")
    print(f"  节点1位移: u1={de[0]:.6e} m, v1={de[1]:.6e} m, w1={de[2]:.6e} m")
    print(f"  节点2位移: u2={de[3]:.6e} m, v2={de[4]:.6e} m, w2={de[5]:.6e} m")
    print(f"  位移列阵 de = {de}")
    
    # ========== 计算刚度矩阵 ==========
    try:
        Ke, L_orig, (cx, cy, cz) = truss3d_element_stiffness(x1, x2, E, A)
    except ValueError as e:
        print(f"\n错误: {e}")
        exit(1)
    
    # 输出几何信息
    print(f"\n【单元几何 (变形前)】")
    print(f"  单元长度 L = {L_orig:.6f} m")
    print(f"  方向余弦 (cx, cy, cz) = ({cx:.6f}, {cy:.6f}, {cz:.6f})")
    
    # 输出刚度矩阵
    print(f"\n【单元刚度矩阵 Ke (6x6)】")
    np.set_printoptions(precision=4, suppress=True, linewidth=100)
    print(Ke)
    
    # ========== 检查刚度矩阵性质 ==========
    print(f"\n【刚度矩阵性质检查】")
    checks = check_stiffness_matrix(Ke)
    
    print(f"  对称性: {'✓ 通过' if checks['对称性'] else '✗ 不通过'}")
    print(f"  对角元素: {'✓ 全部非负' if checks['对角元素非负'] else '✗ 存在负值'}")
    print(f"    对角元素值: {checks['对角元素值']}")
    print(f"  特征值: {checks['特征值']}")
    print(f"  最小特征值: {checks['最小特征值']:.4e}")
    print(f"  零特征值个数: {checks['零特征值个数']}")
    print(f"  矩阵秩: {checks['秩']} (完整秩应为6)")
    print(f"  秩亏损: {checks['秩亏损']}")
    print(f"  行列式: {checks['行列式']:.4e}")
    print(f"  是否奇异: {'是 (存在刚体位移模式)' if checks['是否奇异'] else '否'}")
    
    # ========== 计算应力应变 (方法1: 通过位移) ==========
    print(f"\n【内力与变形 (方法1: 基于位移)】")
    
    # 计算应变 (通过位移)
    delta_L_displacement = (de[3] - de[0]) * cx + (de[4] - de[1]) * cy + (de[5] - de[2]) * cz
    epsilon_displacement = delta_L_displacement / L_orig
    sigma_displacement = E * epsilon_displacement
    N_displacement = sigma_displacement * A
    
    print(f"  轴向应变 ε = {epsilon_displacement:.6e}")
    print(f"  轴向应力 σ = {sigma_displacement/1e6:.4f} MPa")
    print(f"  轴力 N = {N_displacement:.2f} N")
    
    # ========== 计算应力应变 (方法2: 直接通过变形前后长度) ==========
    print(f"\n【内力与变形 (方法2: 基于变形前后长度)】")
    epsilon_direct, sigma_direct, N_direct, L_orig2, L_def, dir_cos = \
        compute_stress_force_from_coords(x1, x2, x1_prime, x2_prime, E, A)
    
    print(f"  变形前长度 L0 = {L_orig2:.6f} m")
    print(f"  变形后长度 L = {L_def:.6f} m")
    print(f"  轴向应变 ε = {epsilon_direct:.6e}")
    print(f"  轴向应力 σ = {sigma_direct/1e6:.4f} MPa")
    print(f"  轴力 N = {N_direct:.2f} N")
    
    # ========== 平衡验证 ==========
    Fe = Ke @ de
    print(f"\n【平衡验证】 Fe = Ke * de = {Fe}")
    
    # ========== 输出汇总 ==========
    print("\n" + "=" * 80)
    print("输入数据汇总")
    print("=" * 80)
    print(f"变形前节点1: ({x1[0]}, {x1[1]}, {x1[2]})")
    print(f"变形前节点2: ({x2[0]}, {x2[1]}, {x2[2]})")
    print(f"变形后节点1: ({x1_prime[0]}, {x1_prime[1]}, {x1_prime[2]})")
    print(f"变形后节点2: ({x2_prime[0]}, {x2_prime[1]}, {x2_prime[2]})")
    print(f"弹性模量 E = {E/1e9:.2f} GPa")
    print(f"截面积 A = {A:.4e} m²")
    print("=" * 80)