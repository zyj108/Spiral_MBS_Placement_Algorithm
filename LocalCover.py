import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from scipy.linalg import solve
from matplotlib.patches import Circle


# 每个节点输入为(x, y)
def distance(A, B):
    return np.sqrt(np.square(A[0] - B[0]) + np.square(A[1] - B[1]))


# 两点中间点位置
def cal_center(p0, p1):
    ux = (p0[0] + p1[0]) / 2
    uy = (p0[1] + p1[1]) / 2
    center = (ux, uy)
    return center

# P_sec最好快排后添加
def one_center_sub(r, P_prio_raw, P_sec_raw):    ##################################如果打算每次只传入一个点，那么P_sec也要写成[point]的形式
    P_prio = P_prio_raw.copy()
    P_sec = P_sec_raw.copy()
    assert P_prio is not P_prio_raw
    assert P_sec is not P_sec_raw

    def distance(A, B):
        return np.sqrt(np.square(A[0] - B[0]) + np.square(A[1] - B[1]))

    # 两点中间点位置
    def cal_center(p0, p1):
        ux = (p0[0] + p1[0]) / 2
        uy = (p0[1] + p1[1]) / 2
        center = (ux, uy)
        return center

    def choose_two_points_and_center(r, points):
        if len(points) >= 2:
            max_dist = 0
            for i, p0_ in enumerate(points):
                for j, p1_ in enumerate(points):
                    if j>i:
                        dist_p0_p1 = distance(p0_, p1_)
                        if dist_p0_p1 >= max_dist:
                            max_dist = dist_p0_p1
                            p0 = p0_
                            p1 = p1_
            center = cal_center(p0, p1)
            # assert distance(p0, p1) < 2 * r
        else:  # 点数如果小于2，就只有一个点，以该点为圆心
            p0 = points[0]
            p1 = points[0]
            center = points[0]
        return p0, p1, center

    def points_outside_circle(center, r_d0_d1, P_sec):
        points_outside = []
        for point in P_sec:
            if distance(center, point) > r_d0_d1:
                points_outside.append(point)
        return points_outside

    def count_degree(p0, p1, p2):
        dist0 = distance(p1, p2)
        dist1 = distance(p0, p2)
        dist2 = distance(p0, p1)
        degree_p0 = math.degrees(math.acos((dist1 * dist1 + dist2 * dist2 - dist0 * dist0) / (2 * dist1 * dist2)))
        degree_p1 = math.degrees(math.acos((dist0 * dist0 + dist2 * dist2 - dist1 * dist1) / (2 * dist0 * dist2)))
        degree_p2 = math.degrees(math.acos((dist1 * dist1 + dist0 * dist0 - dist2 * dist2) / (2 * dist1 * dist0)))

        return [degree_p0, degree_p1, degree_p2]

    def get_center_r_of_circle(A, B, C):
        xa, ya = A[0], A[1]
        xb, yb = B[0], B[1]
        xc, yc = C[0], C[1]
        # 两条边的中点
        x1, y1 = (xa + xb) / 2.0, (ya + yb) / 2.0
        x2, y2 = (xb + xc) / 2.0, (yb + yc) / 2.0
        # 两条线的斜率
        ka = (yb - ya) / (xb - xa) if xb != xa else None
        kb = (yc - yb) / (xc - xb) if xc != xb else None
        alpha = np.arctan(ka) if ka != None else np.pi / 2
        beta = np.arctan(kb) if kb != None else np.pi / 2
        # 两条垂直平分线的斜率
        k1 = np.tan(alpha + np.pi / 2)
        k2 = np.tan(beta + np.pi / 2)
        # 圆心
        y, x = solve([[1.0, -k1], [1.0, -k2]], [y1 - k1 * x1, y2 - k2 * x2])
        # 半径
        r1 = np.sqrt((x - xa) ** 2 + (y - ya) ** 2)
        return x, y, r1

    def find_C(u, A, D, step_four_point):
        if u[0] - A[0] == 0:
            k = 0
        else:
            k = (u[1] - A[1]) / (u[0] - A[0])
        b = u[1] - k * u[0]
        positive_or_neg = (k * D[0] + b) - D[1]
        if positive_or_neg * ((k * step_four_point[0][0] + b) - step_four_point[0][1]) <= 0:
            return step_four_point[0]
        else:
            return step_four_point[1]
        # 绘图查看
        # x = [0, 100, 0.1]
        # plt.plot(x, [k * xi + b for xi in x])
        # plt.ylim(0, 100)
        # plt.scatter(D[0], D[1])
        # positive_or_neg = (k * D[0] + b) - D[1]
        # positive_or_neg * ((k * P_sec[0][0] + b) - P_sec[0][1]) <= 0
        # plt.scatter(u[0], u[1])
        # plt.scatter(A[0], A[1])
        # plt.scatter(step_four_point[0][0], step_four_point[0][1])
        # plt.scatter(step_four_point[1][0], step_four_point[1][1])


    p0, p1, u = choose_two_points_and_center(r, P_prio)
    #对任选两个点的情况中出现只有两个点相同的情况进行处理
    i = 0
    len_P_sec = len(P_sec)
    if p0 == p1:
        while p0 == p1:
            if i == len_P_sec:
                return 0, p0       #如果所有P_sec都遍历完了，都不能组成两个点，则返回唯一的p0
            p1 = P_sec[i]
            i += 1
            if distance(p0, p1)/2 > r:
                p1 = p0
        #新选取的p1需要加入到P_prio中,并重新计算中心点
        P_sec.pop(P_sec.index(p1))
        P_prio.append(p1)
        u = cal_center(p0, p1)


    #定义 oncenter论文中的234步
    def ont_center_234_step(p0, p1, r, P_prio, P_candidate):
        '''如果该函数可以返回是否可以添加P_candidate的点，返回True或False'''
        u = cal_center(p0, p1)
        r_p0_p1 = distance(p0, p1) / 2
        if r_p0_p1 > r:
            print("p0, p1长度出错")

        # #如果P_candidate != 1可以考虑下面的解法
        # #加入一个禁忌列表：
        # #作用: 如果第一次禁忌列表被添加满（ == 圈外点数）, 清空，
        # #     第二次被添加满，强行退出while循环
        # fill_taboo_count = 0
        # taboo_list = []

        p2 = 0      #如果p2为0，说明运算说成中只有2个点； 如果p2是一个点，说明while中可以进入论文中的step
        for just_once in range(2):
            #便于改成while， 为2是因为运行完第一次后，step3需要返回step2进行更新
            if not p2:
                u = cal_center(p0, p1)
                points_ouside = points_outside_circle(u, r_p0_p1, P_candidate)
                if not points_ouside:   #如果为空就跳出循环
                    break
                p2 = P_candidate.pop()
                # p2 = points_ouside[0]
            # step3
            if p2:  # p2如果为一个点，则可以进入step3
                try:
                    degrees = count_degree(p0, p1, p2)
                except:
                    ppp = [p0[0], p1[0], p2[0]]
                    ppp = np.argsort(ppp)
                    degrees = []
                    for i in ppp:
                        if i == 1:
                            degrees.append(180)
                        else:
                            degrees.append(0)
                    print(degrees)
                step_three_points = [p0, p1, p2]
                for i, degree in enumerate(degrees):
                    if degree>=90:
                        p2_ = step_three_points.pop(i)
                        p0_ = step_three_points[0]
                        p1_ = step_three_points[1]
                        if distance(p0_, p1_)/2 > r:    #如果加入p2这个点使得r‘ > r ，令p2 = 0
                            p2 = 0
                            continue #直接跳出这层for 循环，由if p2 == 0进入下一次while循环,避免对p0, p1的改变
                        p0 = p0_
                        p1 = p1_
                        p2 = 0 #因为需要删除一个点，故令p2 = 0
                if p2 == 0:
                    continue

            # step4
            if p2:  # 确保是三个点进入step4
                ux_, uy_, ur = get_center_r_of_circle(p0, p1, p2)
                if ur > r:  #如果半径越界则退出
                    p2=0
                    continue

                in_circle = 1
                for i in P_prio:
                    if distance((ux_, uy_), i) > r:
                        in_circle = 0
                if in_circle == 0: continue

                u = (ux_, uy_)
                points_ouside_step_4 = points_outside_circle((ux_, uy_), ur, P_candidate)
                if not points_ouside_step_4 or points_ouside[0] in points_ouside_step_4:
                    break
                D = points_ouside_step_4.pop()
                step_four_points = [p0, p1, p2]
                step_four_distance = [distance(D, i) for i in step_four_points]
                max_dist_index = step_four_distance.index(np.array(step_four_distance).max())
                A = step_four_points.pop(max_dist_index)
                C = find_C(u, A, D, step_four_points)
                ux_, uy_, ur = get_center_r_of_circle(A, C, D)
                if ur > r:
                    continue
                p0 = A
                p1 = C
                p2 = D

        #如果P_candidate含多个点，则需要修改下面的规则


        for prio in P_prio:
            if distance(u, prio) > r:
                return [0, u]
        return [1, u]

    return ont_center_234_step(p0, p1, r, P_prio, P_sec)


def LocalCover(wk, r, P_prio_raw, P_sec_raw):
    # 对数据进行备份，避免影响输入进来的P_prio_actural和P_sec_actural
    u = (wk[0], wk[1])
    P_prio = P_prio_raw.copy()
    P_sec = P_sec_raw.copy()
    assert P_prio is not P_prio_raw
    assert P_sec is not P_sec_raw

    P_sec_ = []
    while P_sec: #停止的条件是algorithm2 -2中，P_sec已经没有点与u的距离u小于r
        if set(P_sec_) == set(P_sec):
             break
        P_prio = list(set(P_prio))

        P_sec = list(set(P_sec))
        P_sec_ = P_sec[:]

        tmp_P_sec = []
        for pri in P_prio:
            for sec in P_sec:       #####################################################这一步要用P_prio_raw 吗？
                if distance(pri, sec) <= 2 * r:
                    tmp_P_sec.append(sec)

        P_sec = list(set(P_sec))


        dist_P_sec = []
        tmp_P_sec = []
        for sec in P_sec:       ######################################################这一步要用P_sec_raw 吗？
            dist = distance(u, sec)
            if dist <= r:
                P_prio.append(sec)
            else:
                tmp_P_sec.append(sec)
                dist_P_sec.append(dist)

        #取代快排的作用
        # argsort_dist_P_sec = np.argsort(dist_P_sec)
        # tmp_P_sec = np.array(tmp_P_sec)[argsort_dist_P_sec]
        # P_sec = [(two_point[0],two_point[1]) for two_point in np.array(tmp_P_sec)]
        np.random.shuffle(P_sec)

        P_prio = list(set(P_prio))
        P_sec = list(set(P_sec))


        for sec in P_sec:
            if sec in P_prio:
                continue
            whether_append, u_ = one_center_sub(r, P_prio, [sec])
            if whether_append:
                u = u_
                P_prio.append(sec)
                for i in P_prio:
                    if distance(i, u) > r:
                        P_prio.pop(P_prio.index(i))
                P_prio = list(set(P_prio))

    P = P_prio_raw + P_sec_raw
    tmp_P_prio = []
    for point in P:
        if distance(point, u) <= r:
            tmp_P_prio.append(point)
    return u, tmp_P_prio

def draw_circle(x, y, r, alpha):
    ax = plt.gca()
    cir = Circle(xy=(x, y), radius=r, alpha=alpha)
    ax.add_patch(cir)
    ax.plot()

def one_center(r, Points, find_best=False):
    '''输入： 半径， 所要寻找的点
       返回： 圆的中心点，所包含的点
       如果要寻找最优，则find_best=True'''
    if not find_best:
        P_prio_raw = [Points[np.random.choice(range(len(Points)))]]
        wk = P_prio_raw[0]
        u, P_prio = LocalCover(wk, r, P_prio_raw, Points)
    else:
        P_prio = []
        for i in range(len(Points)):
            P_prio_raw = [Points[i]]
            wk = P_prio_raw[0]
            u_, P_prio_ = LocalCover(wk, r, P_prio_raw, Points)
            if len(P_prio_) >= len(P_prio):
                P_prio = P_prio_[:]
                u = u_
    return u, P_prio

if __name__ == '__main__':
    for i in range(1):
        # df = pd.read_csv('./2019-05-29_16_00_00.csv')

        # 生成临时点，用于测试
        points = np.random.randint(0, 100, 50)
        print(points.__repr__())
       #  points = [87,  7, 41, 41, 34, 15, 71, 87, 55,  0, 80, 35, 13, 29, 62, 87, 58,
       # 57, 13, 43, 68, 82, 76, 77, 76, 65, 68, 51, 74, 79, 45, 59, 65, 74,
       # 68, 22, 70, 66, 79, 73, 45, 92,  7, 78, 62, 81, 51, 39, 83, 51]
        data = np.mat(points).reshape(-1, 2)
        df = pd.DataFrame(data=data, columns=['lng', 'lat'])
        lng_lat = [i for i in zip(df['lng'], df['lat'])]
        np.random.randn()


        r = 15
        u = (50, 50)  # 新的MBS要放置的位置
        P_prio_raw = [(35, 40), (40, 30), (30, 40)]

        for point1 in P_prio_raw:
            for point2 in P_prio_raw:
                dist = np.sqrt(np.square(point1[0] - point2[0]) + np.square(point1[1] - point2[1]))
                if dist > r:
                    print('P_prio_raw 设置有误', point1, point2, dist)

        P_prio = P_prio_raw
        P_sec = lng_lat
        # plt.scatter(u[0], u[1])
        # plt.scatter([i[0] for i in P_prio], [i[1] for i in P_prio])
        # plt.scatter([i[0] for i in P_sec], [i[1] for i in P_sec])
        # plt.ylim(0, 100)
        # plt.xlim(0, 100)
        # plt.show()


        wk = P_prio[0]
        # u, P_prio = LocalCover(wk, r,P_prio, P_sec)
        u, P_prio = one_center(r, P_sec, find_best=True)

        #画圆
        cir = draw_circle(u[0], u[1], r, alpha=0.5)
        plt.scatter([i[0] for i in P_sec], [i[1] for i in P_sec], alpha = 0.5, color='g',  s=80)
        #画出圆心
        # plt.scatter(u[0], u[1])
        #画出最开始的P_prio
        plt.scatter([i[0] for i in P_prio_raw], [i[1] for i in P_prio_raw], alpha = 0.5, color='r', s=10)
        plt.scatter([i[0] for i in P_prio], [i[1] for i in P_prio], alpha = 0.5, color='b', s=250)
        plt.axis("equal")
        plt.xlim(0, 100)
        plt.ylim(0, 100)
        plt.show()
        print
        # plt.savefig('figure{}'.format(i))
        # plt.close()
