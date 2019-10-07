
from matplotlib import  pyplot as plt
import numpy as np
import pandas as pd
import copy
import LocalCover
import math
import UAVByKMean
class SpiralMBSPlacement:

    # 通过计算三角形p1p2p3的面积（点在直线左边结果为正，直线右边结果为负）来判断 p3相对于直线p1p2的位置
    def calTri(selt,p1, p2, p3):
        size = p1[0] * p2[1] + p2[0] * p3[1] + p3[0] * p1[1] - p3[0] * p2[1] - p2[0] * p1[1] - p1[0] * p3[1]
        return size

    # 找出据直线最远的点（该点与直线围成的三角形的面积为正且最大）
    def maxSize(self,seq, dot1, dot2, dotSet):
        maxSize = float('-inf')
        maxDot = ()
        online = []
        maxSet = []
        for u in seq:
            size = self.calTri(dot1, dot2, u)
            # 判断点u是否能是三角形u dot1 dot2 的面积为正
            if size < 0:
                continue
            elif size == 0:
                online.append(u)
            # 若面积为正，则判断是否是距离直线最远的点
            if size > maxSize:
                if len(maxDot) > 0:
                    maxSet.append(maxDot)
                maxSize = size
                maxDot = u
            else:
                maxSet.append(u)
        dot2IndexBased = dotSet.index(dot2) + 1

        # 没有分割点,同时可能有点落在直线dot1 dot2上
        if not maxDot:
            # 将直线的点进行一个排序后再按一定的顺寻插入凸包逆时针表示中
            online.sort()
            # 因为是一个逆时针排序的问题，上包不能以dot1为基准进行一个排序，应该以dot2为基准
            # dot1Index = dotSet.index(dot1)
            # for i in range(len(online)):
            #     dotSet.insert(dot1Index,online[i])
            onlineSize = len(online)
            for i in range(onlineSize):
                dotSet.insert(dot2IndexBased + i, online[onlineSize - i - 1])
        # 有分割点
        else:
            # 因为最后凸包是逆时针排序，所以上包部分找到的新的凸包点在数组中都是dot2后面
            dotSet.insert(dot2IndexBased, maxDot)
        return maxSet, maxDot

    # 找出据直线最远的点（该点与直线围成的三角形的面积为负数且最大）
    def minSize(self,seq, dot1, dot2, dotSet):
        minSize = float('inf')
        minDot = ()
        online = []
        minSet = []
        for u in seq:
            size = self.calTri(dot1, dot2, u)
            # 判断点u是否能是三角形u dot1 dot2 的面积为负
            if size > 0:
                continue
            elif size == 0:
                online.append(u)
            # 若面积为负，则判断是否是距离直线最远的点
            if size < minSize:
                if len(minDot) > 0:
                    minSet.append(minDot)
                minDot = u
                minSize = size
            else:
                minSet.append(u)
        # 结果判断
        dot1IndexBased = dotSet.index(dot1) + 1
        # 没找到分割点,同时可能有点落在直线dot1 dot2上
        if not minDot:
            online.sort()
            onlineSize = len(online)
            for i in range(onlineSize):
                dotSet.insert(dot1IndexBased + i, online[i])
        # 有分割点
        else:
            dotSet.insert(dot1IndexBased, minDot)
        return minSet, minDot

    # 上包的递归划分
    def divideUp(self,seq, dot1, dot2, dot3, dot4, dotSet=None):
        # print(dot1, dot2, dot3, dot4)
        # 初始化第一次运行时的参数
        if len(seq) == 0:
            return

        leftSet, rightSet = [], []
        # 划分上包左边的点集
        leftSet, maxDot = self.maxSize(seq, dot1, dot2, dotSet)

        # 对上包左包的点集进一步划分
        if leftSet:
            self.divideUp(leftSet, dot1, maxDot, maxDot, dot2, dotSet)

        # 划分上包右边的点集
        rightSet, maxDot = self.maxSize(seq, dot3, dot4, dotSet)

        # 对上包右包的点集进一步划分
        if rightSet:
            self.divideUp(rightSet, dot3, maxDot, maxDot, dot4, dotSet)
        # return dotSet

    # 下包的递归划分
    def divideDown(self, seq, dot1, dot2, dot3, dot4, dotSet=None):
        # 初始化第一次运行时的参数
        if len(seq) == 0:
            return
        leftSet, rightSet = [], []
        # 划分下包左边的点集
        leftSet, minDot = self.minSize(seq, dot1, dot2, dotSet)

        # 对下包的左包进行进一步划分
        if leftSet:
            self.divideDown(leftSet, dot1, minDot, minDot, dot2, dotSet)

        # 划分下包右包的点集
        rightSet, minDot = self.minSize(seq, dot3, dot4, dotSet)

        # 对下包的右包进一步划分
        if rightSet:
            self.divideDown(rightSet, dot3, minDot, minDot, dot4, dotSet)

    # 递归主函数
    def mainDivide(self,seq):
        # 将序列中的点按横坐标升序排序
        seq.sort()
        res = []
        # 获取横坐标最大、最小的点及横坐标
        dot1 = seq[0]
        dot2 = seq[-1]
        seq1 = []
        maxSize = float('-inf')
        maxDot = ()
        seq2 = []
        minSize = float('inf')
        minDot = ()
        # 对序列划分为直线dot1 dot2左右两侧的点集并找出两个点集的距直线最远点
        for u in seq[1:-1]:
            size = self.calTri(dot1, dot2, u)
            if size > 0:
                if size > maxSize:
                    if len(maxDot) > 0:
                        seq1.append(maxDot)
                    maxSize = size
                    maxDot = u
                    continue
                else:
                    seq1.append(u)
            elif size < 0:
                if size < minSize:
                    if len(minDot) > 0:
                        seq2.append(minDot)
                    minSize = size
                    minDot = u
                    continue
                else:
                    seq2.append(u)

        # 调用内建递归函数
        dotSet = [dot1, dot2]
        if minDot:
            dotSet.insert(1, minDot)
        if maxDot:
            dotSet.append(maxDot)
        self.divideUp(seq1, dot1, maxDot, maxDot, dot2, dotSet)
        self.divideDown(seq2, dot1, minDot, minDot, dot2, dotSet)
        return dotSet

    def execute(self,seq,radius):
        m = 0
        UAVpositionSet = []
        count = 0

        while (seq):
            # print("m=",m)
            # 求凸包并按逆时针顺序排序
            boundarySetInOrder = s.mainDivide(seq)
            innerSet = list(set(seq).difference(set(boundarySetInOrder)))
            # plt.title("overview")
            # plt.scatter([dot[0] for dot in innerSet], [dot[1] for dot in innerSet], color='black')
            # plt.scatter([dot[0] for dot in boundarySetInOrder], [dot[1] for dot in boundarySetInOrder], color='red')
            # x = [dot[0] for dot in res]
            # x.append(res[0][0])
            # y = [dot[1] for dot in res]
            # y.append(res[0][1])
            # plt.plot(x, y)
            # plt.show()
            # plt.plot([dot[0] for dot in boundarySetInOrder], [dot[1] for dot in boundarySetInOrder])
            # plt.show()
            # 在程序刚运行或是上一次循环迭代求得的凸包最后所有的点都被覆盖，
            # 则在当次迭代新求得的凸包中重新随机选择一点，（这里将第一个点作为随机点）
            if  m ==0:
                UAVposition = boundarySetInOrder[0]

            # plt.plot([UAVposition[0]],[UAVposition[1]],"X")
            # plt.show()

            # 求当前UAV的位置，以便在进入下一个迭代的时候按逆时针方向选择其第一个最近的未被覆盖的点
            currentUAVPostion = boundarySetInOrder.index(UAVposition)
            boundarySetInOrderSize = len(boundarySetInOrder)
            nextUAVInd = (currentUAVPostion+1)%boundarySetInOrderSize

            # 第一次调用localCover，优先点为当前无人机的位置（选择的边界点）
            # print(list(set(boundarySetInOrder).difference(set(UAVposition))))
            # boundaryExpectCurrentUAV = list(set(boundarySetInOrder).difference(set(UAVposition)))
            # if count == 3:
            #     print("传入的中心值:",UAVposition,"传入的second:",* list(set(boundarySetInOrder).difference(set(UAVposition))))

            UAVposition,prioSet = LocalCover.LocalCover(UAVposition,radius,[UAVposition],
                                                  list(set(boundarySetInOrder).difference(set(UAVposition))))

            prioBoCopy = copy.deepcopy(prioSet)

            # 第二次调用localCover,优先点为上一次调用的优先结果点集
            # 返回结果即为该次放置无人机的位置以及覆盖的点
            UAVposition,prioSet = LocalCover.LocalCover(UAVposition,radius,prioSet,innerSet)
            # print("中心值：",UAVposition,"优先点：",*prioSet,sep="\t")
            # print("prioSet",*prioSet,sep="\t")
            UAVpositionSet.append(UAVposition)
            count = count+1
            # plt.scatter([dot[0] for dot in innerSet], [dot[1] for dot in innerSet], color='black')
            # plt.scatter([dot[0] for dot in boundarySetInOrder], [dot[1] for dot in boundarySetInOrder], color='red')
            # plt.plot([UAVposition[0]], [UAVposition[1]], "X")
            # theta = np.arange(0, 2 * np.pi, 0.01)
            # x = UAVposition[0] + radius * np.cos(theta)
            # y = UAVposition[1] + radius * np.sin(theta)
            # plt.plot(x, y)
            # plt.show()
            # 计算下一个迭代距离当前选择的UAV位置第一近的未被覆盖的点
            flag = True
            while nextUAVInd != currentUAVPostion:
                if boundarySetInOrder[nextUAVInd] not in prioBoCopy:
                    UAVposition = boundarySetInOrder[nextUAVInd]
                    m = m + 1
                    flag = False
                    break;
                else:
                    nextUAVInd = (nextUAVInd + 1) % boundarySetInOrderSize
            if flag:
                m = 0

            # 去除已覆盖的点，进入下一个无人机位置的计算
            seq  = list(set(seq).difference(set(prioSet)))

        return UAVpositionSet

if __name__ == "__main__":
    s = SpiralMBSPlacement()
    radius = 100
    # seq0 = [(random.randint(0, 50), random.randint(0, 50)) for x in range(50)]
    # seq0 = list(set(seq0))
    # s.execute(seq0)
    data = pd.read_csv("深圳市-欢乐谷.csv",header=0,usecols=[2,3,4]).values
    # data = [(int(dot[0]),int(dot[1])) for dot in data if dot[2]  ]
    data = [(math.modf(dot[0])[1],math.modf(dot[1])[1]) for dot in data   ]
    # data = [(dot[0],dot[1]) for dot in data  ]
    # data = list(set(data))
    # for dot in data:
    #     print(dot[0],dot[1])
    data = list(set(data[0:400]))

    #  使用K-means
    # UAVpositionSet = UAVByKMean.binarysearch(data,  radius)
    # print(*data,sep='\n')
    UAVpositionSet = s.execute(data,radius)
    print(len(UAVpositionSet))
    x = [dot[0]  for dot in data]
    y = [dot[1] for dot in data]
    plt.scatter(x, y, color='black')
    # for a, b in zip(x, y):
    #     plt.text(a, b, (a, b), ha='center', va='bottom', fontsize=10)
    ## 按顺寻连接无人机的位置
    for UAV in UAVpositionSet:
        # print(UAV)
        print("[",UAV[0],",",UAV[1],"]",sep="")
        theta = np.arange(0, 2 * np.pi, 0.01)
        x = UAV[0] + radius * np.cos(theta)
        y = UAV[1] + radius * np.sin(theta)
        plt.plot(x, y)  
    #
    x = [dot[0] for dot in UAVpositionSet]
    y = [dot[1] for dot in UAVpositionSet]
    # # plt.xlim(126880000000,126880000700)
    # # plt.ylin
    # plt.scatter(x, y, color='red')
    plt.plot(x,y,color = "red")
    plt.show()