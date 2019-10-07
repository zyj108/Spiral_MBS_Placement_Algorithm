from sklearn.cluster import KMeans
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import math
def calculate_distance(point1,point2):
    d_x = point1[0] - point2[0]
    d_y =  point1[1] - point2[1]
    #计算两点之间的距离
    distance = math.sqrt(d_x**2 + d_y**2)
    return distance
def binarysearch(data,radius):

    n_clusters = list(range(2,len(data)+1))
    length = len(n_clusters)
    UAVPositionSet = []
    count = 0
    minSize = float('inf')
    aveUAV = 0
    while count < 5 :
        low = 0
        high = length - 1
        UAVPosition = []
        while low <= high:
            mid = int(low + ((high - low) / 2))  ##使用(low+high)/2会有整数溢出的问题
            kmeans = KMeans(n_clusters=n_clusters[mid], max_iter=1000)
            predicted = kmeans.fit_predict(data)
            # 记录簇结果
            cluster_result = {}
            for i in range(n_clusters[mid]):
                cluster_result[i] = list()
            for i in range(len(predicted)):
                cluster_result[predicted[i]].append(data[i])

            centroids = kmeans.cluster_centers_
            outloopFlag = True
            for i in range(n_clusters[mid]):
                innerloopFlag = False
                for point in cluster_result[i]:
                    if calculate_distance(point,centroids[i]) > radius:
                        innerloopFlag = True
                        break
                if innerloopFlag:
                    low = mid+1
                    outloopFlag = False
                    break

            if outloopFlag:
                high = mid - 1
                UAVPosition = centroids

        print("当前UAV的数量；",len(UAVPosition))
        count = count+1
        aveUAV = aveUAV + len(UAVPosition)
        if len(UAVPosition) < minSize:
            minSize = len(UAVPosition)
            UAVPositionSet = UAVPosition

    return UAVPositionSet

