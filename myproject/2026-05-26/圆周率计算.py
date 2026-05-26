import math
import os

# 设置保存路径为桌面
desktop_path = os.path.join(os.path.expanduser("~"), "Desktop")
output_file = os.path.join(desktop_path, "convergence_plot_data.txt")

def method_polygon(n_values):
    """方法一：多边形逼近法"""
    results = []
    exact_pi = math.pi
    
    for n in n_values:
        pi_n = n * math.sin(math.pi / n)
        error = abs(exact_pi - pi_n)
        h = 1.0 / n
        results.append((n, h, error))
    
    return results

def method_richardson_extrapolation(n_values):
    """方法二：Richardson 外推法"""
    pi_values = [n * math.sin(math.pi / n) for n in n_values]
    exact_pi = math.pi
    
    results = []
    
    # 第一点无法外推，使用原始值
    n1 = n_values[0]
    pi1 = pi_values[0]
    error1 = abs(exact_pi - pi1)
    h1 = 1.0 / n1
    results.append((n1, h1, error1))
    
    # 外推后续点
    for i in range(1, len(pi_values)):
        n_prev = n_values[i-1]
        n_curr = n_values[i]
        pi_prev = pi_values[i-1]
        pi_curr = pi_values[i]
        
        r = n_curr / n_prev
        p = 2
        extrap_val = pi_curr + (pi_curr - pi_prev) / (r**p - 1)
        
        error = abs(exact_pi - extrap_val)
        h = 1.0 / n_curr
        results.append((n_curr, h, error))
    
    return results

def save_for_plot(polygon_results, extrapolation_results, filename):
    """保存为绘图数据：h 列，多边形误差列，外推误差列"""
    # 创建外推结果字典
    extrap_dict = {n: (h, error) for n, h, error in extrapolation_results}
    
    with open(filename, 'w', encoding='utf-8') as f:
        f.write('TITLE = "两种方法收敛性对比"\n')
        f.write('VARIABLES = "h", "error_polygon", "error_extrapolation"\n')
        
        # 统计有效数据点数
        valid_points = []
        for n, h, error_poly in polygon_results:
            if n in extrap_dict:
                h_ext, error_ext = extrap_dict[n]
                valid_points.append((h, error_poly, error_ext))
        
        f.write(f'ZONE I={len(valid_points)}, F=POINT\n')
        
        for h, error_poly, error_ext in valid_points:
            f.write(f"{h:15.12e}    {error_poly:15.12e}    {error_ext:15.12e}\n")
    
    print(f"绘图数据已保存: {filename}")

def save_for_tecplot_loglog(polygon_results, extrapolation_results, filename):
    """保存为 log-log 格式（log10(h) vs log10(error)）"""
    extrap_dict = {n: (h, error) for n, h, error in extrapolation_results}
    
    with open(filename, 'w', encoding='utf-8') as f:
        f.write('TITLE = "两种方法收敛率对比 (log-log)"\n')
        f.write('VARIABLES = "log10_h", "log10_error_polygon", "log10_error_extrapolation"\n')
        
        valid_points = []
        for n, h, error_poly in polygon_results:
            if n in extrap_dict:
                h_ext, error_ext = extrap_dict[n]
                if error_poly > 0 and error_ext > 0:
                    log_h = math.log10(h)
                    log_error_poly = math.log10(error_poly)
                    log_error_ext = math.log10(error_ext)
                    valid_points.append((log_h, log_error_poly, log_error_ext))
        
        f.write(f'ZONE I={len(valid_points)}, F=POINT\n')
        
        for log_h, log_error_poly, log_error_ext in valid_points:
            f.write(f"{log_h:12.6f}    {log_error_poly:12.6f}    {log_error_ext:12.6f}\n")

def save_to_readable_table(polygon_results, extrapolation_results, filename):
    """保存为易读表格（只包含 h、多边形误差、外推误差）"""
    extrap_dict = {n: (h, error) for n, h, error in extrapolation_results}
    
    with open(filename, 'w', encoding='utf-8') as f:
        f.write("=" * 70 + "\n")
        f.write("两种方法收敛性对比\n")
        f.write(f"精确 π = {math.pi:.16f}\n")
        f.write("=" * 70 + "\n\n")
        
        f.write(f"{'n':>8} | {'h = 1/n':>15} | {'多边形逼近 e_n':>20} | {'外推法 e_n':>20}\n")
        f.write("-" * 70 + "\n")
        
        for n, h, error_poly in polygon_results:
            if n in extrap_dict:
                h_ext, error_ext = extrap_dict[n]
                f.write(f"{n:8d} | {h:15.12e} | {error_poly:20.12e} | {error_ext:20.12e}\n")
            else:
                f.write(f"{n:8d} | {h:15.12e} | {error_poly:20.12e} | {'':20s}\n")
        
        f.write("=" * 70 + "\n")

def print_console(polygon_results, extrapolation_results):
    """控制台输出"""
    extrap_dict = {n: (h, error) for n, h, error in extrapolation_results}
    
    print("\n" + "=" * 70)
    print("两种方法收敛性对比")
    print(f"精确 π = {math.pi:.16f}")
    print("=" * 70)
    print(f"{'n':>8} | {'h = 1/n':>15} | {'多边形逼近 e_n':>20} | {'外推法 e_n':>20}")
    print("-" * 70)
    
    for n, h, error_poly in polygon_results:
        if n in extrap_dict:
            h_ext, error_ext = extrap_dict[n]
            print(f"{n:8d} | {h:15.12e} | {error_poly:20.12e} | {error_ext:20.12e}")
        else:
            print(f"{n:8d} | {h:15.12e} | {error_poly:20.12e} | {'':20s}")

def main():
    # n 值序列
    n_values = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192]
    
    # 计算两种方法
    polygon_results = method_polygon(n_values)
    extrapolation_results = method_richardson_extrapolation(n_values)
    
    # 控制台输出
    print_console(polygon_results, extrapolation_results)
    
    # 保存到桌面
    plot_file = os.path.join(desktop_path, "convergence_plot_data.dat")
    loglog_file = os.path.join(desktop_path, "convergence_loglog.dat")
    table_file = os.path.join(desktop_path, "convergence_table.txt")
    
    save_for_plot(polygon_results, extrapolation_results, plot_file)
    save_for_tecplot_loglog(polygon_results, extrapolation_results, loglog_file)
    save_to_readable_table(polygon_results, extrapolation_results, table_file)
    
    print("\n" + "=" * 70)
    print(f"文件已保存到桌面：")
    print(f"  1. {plot_file} ← h, error_polygon, error_extrapolation")
    print(f"  2. {loglog_file} ← log10(h), log10(error)")
    print(f"  3. {table_file} ← 易读表格")
    print("=" * 70)

if __name__ == "__main__":
    main()