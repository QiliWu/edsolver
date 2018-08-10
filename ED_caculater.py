# -*- coding: utf-8 -*-
from math import sqrt, sin, cos, acos, pi

def crystal_system(a, b, c, alpha, beta, gamma):
    #根据输入的晶胞参数值判断物相所属晶系
    if alpha == beta == gamma == 90:
        if a == b == c:
            return 'cubic'
        elif a != b and b != c and a != c:
            return 'orthorhombic'
        else:
            return 'tetragonal'
    elif a == b and alpha == beta == 90 and gamma == 120:
        return 'hexagonal'
    elif a != b and b != c and a != c and alpha != beta and beta != gamma and alpha != gamma:
        return 'triclinic'
    elif a != b and b != c and a != c and alpha == gamma == 90 != beta:
        return 'monoclinic'
    else:
        print('请输入有效格式的晶胞参数！')
        return None

def d_caculate(a, b, c, alpha, beta, gamma, h, k, l, crys_sys):
    #计算指定晶面的晶面间距
    if crys_sys:
        try:
            alpha = alpha * pi / 180
            beta = beta * pi / 180
            gamma = gamma * pi / 180
            square_g = 0
            if crys_sys == 'cubic':
                square_g = (h*h + k*k + l*l)/(a*a)
            elif crys_sys == 'orthorhombic':
                square_g = (h/a)**2 + (k/b)**2 + (l/c)**2
            elif crys_sys == 'tetragonal':
                square_g = (h*h + k*k)/(a*a) + (l/c)**2
            elif crys_sys == 'hexagonal':
                square_g = 4 * (h*h + h*k + k*k)/(3*a*a) + (l/c)**2
            elif crys_sys == 'triclinic':
                square_v = ((a*b*c)**2)*(1-cos(alpha)*cos(alpha)-cos(beta)*cos(beta)-
                                         cos(gamma)*cos(gamma) +2*cos(alpha)*cos(beta)*cos(gamma))
                square_g = (1/square_v)*\
                           ((h*b*c*sin(alpha))**2 + (k*a*c*sin(beta))**2 + (l*b*a*sin(gamma))**2
                            + 2*h*k*a*b*(c**2)*(cos(alpha)*cos(beta)-cos(gamma))
                            + 2*k*l*b*c*(a**2)*(cos(gamma)*cos(beta)-cos(alpha))
                            + 2*l*h*a*c*(b**2)*(cos(alpha)*cos(gamma)-cos(beta)))
            elif crys_sys == 'monoclinic':
                square_g = (h/(a*sin(beta)))**2 + (k/b)**2 + (l/(c*sin(beta)))**2 \
                           - 2*h*l*cos(beta)/(a*c*(sin(beta))**2)
            return round(1/sqrt(square_g), 3)
        except:
            print('晶胞参数输入有误，请检查后重新尝试！')
            return None
    else:
        return None

def theta_caculate(a, b, c, alpha, beta, gamma, h_1, k_1, l_1, h_2, k_2, l_2, crys_sys):
    #计算两个晶面间的夹角
    try:
        d_1 = d_caculate(a, b, c, alpha, beta, gamma, h_1, k_1, l_1, crys_sys)
        d_2 = d_caculate(a, b, c, alpha, beta, gamma, h_2, k_2, l_2, crys_sys)
        alpha = alpha * pi / 180
        beta = beta * pi / 180
        gamma = gamma * pi / 180
        cos_theta = 0
        if crys_sys == 'cubic':
            cos_theta = (h_1*h_2 + k_1*k_2 + l_1*l_2) * d_1 * d_2 / (a**2)
        elif crys_sys == 'orthorhombic':
            cos_theta = (h_1*h_2/(a*a) + k_1*k_2/(b*b) + l_1*l_2/(c*c)) * d_1 * d_2
        elif crys_sys == 'tetragonal':
            cos_theta = ((h_1*h_2 + k_1*k_2)/(a*a) + l_1*l_2/(c*c)) * d_1 * d_2
        elif crys_sys == 'hexagonal':
            cos_theta = (4/(3*a*a))*(h_1*h_2 + k_1*k_2 + 0.5*(h_1*k_2+h_2*k_1) + l_1*l_2/(c*c)) * d_1 * d_2
        elif crys_sys == 'triclinic':
            square_v = ((a * b * c) ** 2) * (1 - cos(alpha) * cos(alpha) - cos(beta) * cos(beta)
                                             - cos(gamma) * cos(gamma) + 2 * cos(alpha) * cos(beta) * cos(gamma))
            A = (h_1*h_2*(b*c*sin(alpha))**2 + k_1*k_2*(a*c*sin(beta))**2 + l_1*l_2*(a*b*sin(gamma))**2
                 + a*b*c*c*(cos(alpha)*cos(beta) - cos(gamma))*(k_1*h_2 + h_1*k_2)
                 + a*b*b*c*(cos(alpha)*cos(gamma) - cos(beta))*(h_1*l_2 + h_2*l_1)
                 + a*a*b*c*(cos(beta)*cos(gamma) - cos(alpha))*(k_1*l_2 + k_2*l_1)
                 )/square_v
            cos_theta = A * d_1 * d_2
        elif crys_sys == 'monoclinic':
            cos_theta = (h_1*h_2/((a*sin(beta))**2) + k_1*k_2/(b*b) + (l_1*l_2/((c*sin(beta))**2))
                         - (h_2*l_1 + h_1*l_2)*cos(beta)/(a*c*(sin(beta)**2))) \
                        * d_1 * d_2
        return round(acos(cos_theta)*180/pi, 2)
    except:
        return None

def uvw_caculate(h_1, k_1, l_1, h_2, k_2, l_2):
    #计算两个晶面的晶带轴
    u = k_1*l_2 - k_2*l_1
    v = l_1*h_2 - l_2*h_1
    w = h_1*k_2 - k_1*h_2
    #去除公约数，获得最简化组合
    for i in range(max(map(abs, [u,v,w])), 1, -1):
        if u%i == 0 and v%i ==0 and w%i == 0:
            u, v, w = u/i, v/i, w/i
    return (int(u),int(v), int(w))

def third_spot_info(h_1, k_1, l_1, h_2, k_2, l_2):
    #根据平行四边形法则计算第三个衍射点对应的晶面指数
    h_3 = h_2 - h_1
    k_3 = k_2 - k_1
    l_3 = l_2 - l_1
    return (h_3,k_3,l_3)

def d_to_hkl(a, b, c, alpha, beta, gamma, d, delta_d, crys_sys):
    #根据晶面间距推算在误差范围内符合的晶面
    up = d + delta_d
    low = d - delta_d
    result = []
    h, k, l = 1, 1, 1
    try:
        #依次计算h,k,l的边界值
        while d_caculate(a, b, c, alpha, beta, gamma, h, 0, 0, crys_sys) >= low:
            h = h + 1
        while d_caculate(a, b, c, alpha, beta, gamma, 0, k, 0, crys_sys) >= low:
            k = k + 1
        while d_caculate(a, b, c, alpha, beta, gamma, 0, 0, l, crys_sys) >= low:
            l = l + 1
        for i in range(h-1, -h, -1):
            for j in range(k-1, -k,-1):
                for p in range(l-1, -l, -1):
                    if i == j == p == 0:
                        continue
                    elif d_caculate(a, b, c, alpha, beta, gamma, i, j, p, crys_sys) >= low and \
                                    d_caculate(a, b, c, alpha, beta, gamma, i, j, p, crys_sys) <= up:
                        result.append((i,j,p))
        if result:
            return result
        else:
            print('找不到d值对应晶面，请尝试增加delta_d和delta_theta')
            return None
    except:
        return None

def find_hkl(a, b, c, alpha, beta, gamma, d_1, d_2, delta_d, theta, delta_theta, crys_sys):
    #根据晶胞参数，相邻两个衍射点的晶面间距值，它们与中心点构成的夹角，以及给定的误差范围，计算满足条件的衍射点晶面指数组合
    result = []
    hkl_1_list = d_to_hkl(a, b, c, alpha, beta, gamma, d_1, delta_d, crys_sys)
    if hkl_1_list:
        hkl_2_list = d_to_hkl(a, b, c, alpha, beta, gamma, d_2, delta_d, crys_sys)
        if hkl_2_list:
            for s in hkl_1_list:
                h_1, k_1, l_1 = s[0], s[1], s[2]
                for p in hkl_2_list:
                    h_2, k_2, l_2 = p[0], p[1], p[2]
                    if h_1*k_2==h_2*k_1 and h_1*l_2==h_2*l_1 and k_1*l_2==k_2*l_1:
                        #两个晶面平行，不符合
                        continue
                    elif ((-h_1,-k_1,-l_1), (-h_2, -k_2, -l_2)) in result:
                        #该组合已经在结果里了，只是取的是反向晶面指数
                        continue
                    else:
                        if theta_caculate(a, b, c, alpha, beta, gamma, h_1, k_1, l_1, h_2, k_2, l_2, crys_sys) <= theta + delta_theta \
                                and theta_caculate(a, b, c, alpha, beta, gamma, h_1, k_1, l_1, h_2, k_2, l_2, crys_sys) >= theta - delta_theta:
                            result.append(((h_1,k_1,l_1), (h_2, k_2, l_2)))
            if result:
                return result
            else:
                print("找不到匹配的h1k1l1和h2k2l2组合, 请尝试增加Δd和Δθ")
                return None
        else:
            return None

def type_check(args):
    # 检查用户输入值得数据有效性
    for i in args:
        if type(args[i]) != type(1) and type(args[i]) != type(1.2):
            print('请输入整数或小数格式的{}'.format(i))
            return False
        elif args[i] <= 0:
            print('请输入一个正数{}'.format(i))
            return False
    return True

def angle_check(theta_1, d_1, d_2, delta_d, delta_theta):
    #添加一个角度值范围限定函数
    if theta_1 > 90:
        print('请输入一个锐角的theta_1')
        return False
    elif delta_d/d_1 > 0.3 or delta_d/d_2 > 0.3:
        print('d值允许误差超过30%, 结果不准确！')
        return False
    elif delta_theta/theta_1 > 0.1:
        print('theta值允许误差超过30%, 结果不准确！')
        return False
    return True

def main(a=None, b=None, c=None, alpha=None, beta=None, gamma=None,
           d_1=None, d_2=None, delta_d=None, theta_1=None, theta_2=None, delta_theta=None):
    # 首先判断输入值类型，正确了再进行下一步
    args = {'a': a, 'b': b, 'c': c, 'alpha': alpha, 'beta': beta, 'gamma': gamma, 'd_1': d_1,'d_2': d_2,
         'delta_d': delta_d, 'theta_1': theta_1, 'theta_2': theta_2, 'delta_theta': delta_theta}
    if type_check(args) and angle_check(theta_1, d_1, d_2, delta_d, delta_theta):
        crys_sys = crystal_system(a, b, c, alpha, beta, gamma)
        if crys_sys:
            hkl_ls = find_hkl(a, b, c, alpha, beta, gamma,d_1, d_2, delta_d, theta_1, delta_theta, crys_sys)
            result = {}
            saved_d1 = []
            saved_d2 = []
            saved_d3 = []
            if hkl_ls:
                n = 1
                for i in hkl_ls:
                    h_1, k_1, l_1, h_2, k_2, l_2 = i[0][0], i[0][1], i[0][2], i[1][0], i[1][1], i[1][2]
                    d_1_theoretical = d_caculate(a, b, c, alpha, beta, gamma, h_1, k_1, l_1, crys_sys)
                    d_2_theoretical = d_caculate(a, b, c, alpha, beta, gamma, h_2, k_2, l_2, crys_sys)
                    theta_1_theoretical = theta_caculate(a, b, c, alpha, beta, gamma, h_1, k_1, l_1, h_2, k_2, l_2, crys_sys)
                    delta_theta_1 = round(abs(theta_1_theoretical-theta_1), 2)
                    uvw = uvw_caculate(h_1, k_1, l_1, h_2, k_2, l_2)
                    hkl_3 = third_spot_info(h_1, k_1, l_1, h_2, k_2, l_2)
                    d_3 = d_caculate(a, b, c, alpha, beta, gamma, hkl_3[0], hkl_3[1], hkl_3[2], crys_sys)
                    theta_2_theoretical = theta_caculate(a, b, c, alpha, beta, gamma, h_2, k_2, l_2, hkl_3[0], hkl_3[1], hkl_3[2], crys_sys)

                    if abs(theta_2_theoretical-theta_2) <= delta_theta:
                        if d_1_theoretical in saved_d1 and d_2_theoretical in saved_d2 and d_3 in saved_d3:
                            #去除一些等价的组合
                            continue
                        else:
                            result[n] = {'d1': d_1_theoretical, 'd2': d_2_theoretical, 'd3': d_3,
                                     'h1k1l1': (h_1, k_1, l_1), 'h2k2l2': (h_2, k_2, l_2),'h3k3l3': hkl_3,
                                     'uvw': uvw, 'theta_1': theta_1_theoretical, 'theta_2': theta_2_theoretical,
                                     'delta_theta_1': delta_theta_1}
                            n = n + 1
                            saved_d1.append(d_1_theoretical)
                            saved_d2.append(d_2_theoretical)
                            saved_d3.append(d_3)
                if result:
                    return result
                else:
                    print('没有找到匹配的结果，请尝试增加delta_d或delta_theta')
                    return '抱歉！！！'
            else:
                return 'Sorry!!!'
    else:
        return ''

print(main(a=12.163, b=3.735, c=6.513, alpha=90, beta=127.29, gamma=90,
           d_1=3.281, d_2=3.091, delta_d=0.3, theta_1=62.27, theta_2=58.05, delta_theta=2))

