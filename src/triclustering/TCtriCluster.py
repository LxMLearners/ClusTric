########################################################################
#                                                                      #
#               TCtriCluster : *Python Version*                        #
#                                                                      #
#                                                                      #
#               Author: Diogo F Soares                                 #
#                                                                      #
#        Adapted from       triCluster (Zhao and Zaki)                 #
#                                                                      #
########################################################################

from collections import namedtuple
import copy
import itertools
from sortedcontainers import SortedSet

Name = namedtuple("Name", "name")


def fZero(f: float):
    return abs(f) < 0.0000001


BRIEF_OUT = False


class Ratio:
    def __init__(self):
        self.val = None
        self.trans = None
        self.f1Neg = None

    def __lt__(self, other):
        return self.val < other.val


class Array3D:

    def __init__(self, t: int, s: int, g: int):
        self.T = t
        self.S = s
        self.G = g
        self._data = dict()
        self.tName = [None] * self.T  # list of Name()
        self.sName = [None] * self.S
        self.gName = [None] * self.G

    def dat(self, t: int, s: int, g: int):
        return self._data[(t, s, g)]
        # return self._data[t*(S*G)+s*G+g] ???

    def setdat(self, t: int, s: int, g: int, val):
        self._data[(t, s, g)] = val

    def clear(self):
        del self._data
        del self.tName
        del self.sName
        del self.gName

    def show(self):
        print("Total Times:\t", self.T)
        print("Total Samples:\t", self.S)
        print("Total Genes:\t", self.G)

        i = 0
        while i < self.T:
            print("Time\t", i)
            print("ID\tNAME")
            k = 0
            while k < self.S:
                print("\tS-", k)
                k += 1

            j = 0
            while j < self.G:
                print(j, "\tG-", j)
                k = 0
                while k < self.S:
                    print("\t{}".format(self.dat(i, k, j)))
                    k += 1
                j += 1
            i += 1


class TabInput:
    def open(self, file):
        self._filePointer = open(file, 'r')

    def close(self):
        self._filePointer.close()

    def getLine(self):
        self._line = self._filePointer.readline()
        self._sub = ""
        self.k = 0

    def get2tab(self):
        p = 0
        if self._line[self.k] == '\t' or self._line[self.k] == '\n':
            self._sub += ''
            if self._line[self.k] == '\t':
                self.k += 1
        else:
            while True:
                self._sub += self._line[self.k]
                #p += 1
                self.k += 1
                if not (self._line[self.k] != '\t' and self._line[self.k] != '\n' and self.k < len(self._line)):
                    break
            self.k += 1
            self._sub += ''

        s = self._sub
        self._sub = ""
        return s


class Input:

    def __init__(self, filename):
        self._file = filename
        self._in = TabInput()

    def readTabFile(self):

        self._in.open(self._file)

        self._in.getLine()
        self._in.get2tab()
        p = self._in.get2tab()  # T
        T = int(p)

        self._in.getLine()
        self._in.get2tab()
        p = self._in.get2tab()  # S
        S = int(p)

        self._in.getLine()
        self._in.get2tab()
        p = self._in.get2tab()  # G
        G = int(p)

        self._pArray3d = Array3D(T, S, G)

        i = 0
        while i < T:
            self._in.getLine()
            self._in._line
            self._in.get2tab()
            p = self._in.get2tab()
            self._pArray3d.tName[i] = p
            self._in.getLine()

            if i == 0:
                self._in.get2tab()
                self._in.get2tab()  # ID AND NAME

                k = 0
                while k < S:
                    p = self._in.get2tab()
                    self._pArray3d.sName[k] = p
                    k += 1

            j = 0
            while j < G:
                self._in.getLine()
                self._in.get2tab()
                p = self._in.get2tab()  # ID AND NAME
                if i == 0:
                    self._pArray3d.gName[j] = p

                k = 0
                while k < S:
                    p = self._in.get2tab()
                    self._pArray3d.setdat(i, k, j, float(''.join(itertools.takewhile(
                        lambda x: str.isdigit(x) or x == '-' or x == '.', p))))
                    k += 1
                j += 1
            i += 1

        return self._pArray3d


class Rect:
    def __init__(self):
        self.deleted = False
        self.S = SortedSet()
        self.G = SortedSet()

    def clear(self):
        self.G.clear()
        self.S.clear()

    def contain(self, r):
        """
            r is type of Rect
        """
        if len(self.S) < len(r.S) or len(self.G) < len(r.G):
            return False

        U = self.S.union(r.S)
        if len(U) > len(self.S):
            return False

        U = self.G.union(r.G)
        if len(U) > len(self.G):
            return False

        return True

    def copy(self, r):
        #r = Rect()
        self.S = copy.deepcopy(r.S)
        self.G = copy.deepcopy(r.G)
        self.deleted = copy.deepcopy(r.deleted)

    def containPoint(self, s: int, g: int):
        if s not in self.S:
            return False
        if g not in self.G:
            return False
        return True

    def elementNum(self):
        return len(self.S) * len(self.G)

    def show(self, array3d: Array3D, T: int):
        sFinal = "|S| x |G|: {} x {} \n".format(len(self.S), len(self.G))
        sFinal += "                            "
        for its in self.S:
            sFinal += '{:>14}'.format(str(array3d.sName[its]))
        sFinal += '\n'
        for itg in self.G:
            sFinal += '{:>30}\t'.format(str(array3d.gName[itg]))
            for its in self.S:
                sFinal += "      {:6.2f}  ".format(array3d.dat(T, its, itg))
            sFinal += '\n'

        return sFinal


class Rects:

    def __init__(self):
        self.rectVec = list()

    def clear(self):
        for rec in self.rectVec:
            rec.clear()
        self.rectVec.clear()

    def contain(self, r: Rect):
        for rec in self.rectVec:
            if rec.contain(r):
                return True
        return False

    def del_becontained(self, r: Rect):
        for rec in self.rectVec:
            if r.contain(rec):
                rec.S.clear()
                rec.G.clear()
                self.rectVec.remove(rec)

    def addIfMax(self, r: Rect):
        if not self.contain(r):
            self.del_becontained(r)
            self.rectVec.append(r)

    def validNum(self):
        num = 0
        for rec in self.rectVec:
            if not rec.deleted:
                num += 1
        return num

    def show(self, array3d: Array3D, T: int):
        i = 0
        print("================================================ Time: ",
              T, " =================================================")

        for it in self.rectVec:
            if not it.deleted:
                print("Rect ", i, ": ")
                i += 1
                print(it.show(array3d, T))


class Cube:

    def __init__(self):
        self.T = SortedSet()
        self.S = SortedSet()
        self.G = SortedSet()
        self.deleted = False

    def __eq__(self, other):
        if other is None:
            return False

        if not isinstance(other, Cube):
            return False

        return self.T == other.T and self.S == other.S and self.G == other.G and self.deleted == other.deleted

    def clear(self):
        self.G.clear()
        self.S.clear()
        self.T.clear()

    def copy(self, c):
        self.T = copy.deepcopy(c.T)
        self.S = copy.deepcopy(c.S)
        self.G = copy.deepcopy(c.G)
        self.deleted = False

    def contain(self, c):
        """
            r is type of Rect
        """
        if len(self.T) < len(c.T) or len(self.S) < len(c.S) or len(self.G) < len(c.G):
            return False

        U = self.T.union(c.T)
        if len(U) > len(self.T):
            return False

        U = self.S.union(c.S)
        if len(U) > len(self.S):
            return False

        U = self.G.union(c.G)
        if len(U) > len(self.G):
            return False

        return True

    def containPoint(self, t: int, s: int, g: int):
        if t not in self.T:
            return False
        if s not in self.S:
            return False
        if g not in self.G:
            return False
        return True

    def elementNum(self):
        return len(self.T) * len(self.S) * len(self.G)

    def show(self, array3d: Array3D):
        sFinal = ""

        sFinal += "\n|T|x|S|x|G|: " + \
            str(len(self.T)) + "x" + str(len(self.S)) + \
            "x" + str(len(self.G)) + '\n'

        for itt in self.T:
            sFinal += "    Time: {}\n".format(array3d.tName[itt])
            sFinal += "                            "
            for its in self.S:
                sFinal += '{:>14}'.format(str(array3d.sName[its]))
            sFinal += '\n'
            for itg in self.G:
                sFinal += '{:>30}\t'.format(str(array3d.gName[itg]))
                for its in self.S:
                    sFinal += "      {:6.2f}  ".format(
                        array3d.dat(itt, its, itg))
                sFinal += '\n'
            sFinal += '\n'

        return sFinal

    def variance(self, dimension, array3d):
        ave = 0.0
        sum = 0.0
        SUM = 0.0
        if dimension == 'T':
            it1 = self.S
            it2 = self.G
            it3 = self.T
        elif dimension == 'S':
            it1 = self.T
            it2 = self.G
            it3 = self.S
        else:
            it1 = self.T
            it2 = self.S
            it3 = self.G

        for i1 in it1:
            for i2 in it2:
                ave = 0.0
                for i3 in it3:
                    if dimension == 'T':
                        ave += round(array3d.dat(i3, i1, i2), 2)
                    elif dimension == 'S':
                        ave += round(array3d.dat(i1, i3, i2), 2)
                    else:
                        ave += round(array3d.dat(i1, i2, i3), 2)

                if dimension == 'T':
                    ave /= len(self.T)
                elif dimension == 'S':
                    ave /= len(self.S)
                else:
                    ave /= len(self.G)

                sum = 0.0
                for i3 in it3:
                    if dimension == 'T':
                        sum += abs(array3d.dat(i3, i1, i2) - ave)
                    elif dimension == 'S':
                        sum += abs(array3d.dat(i1, i3, i2) - ave)
                    else:
                        sum += abs(array3d.dat(i1, i2, i3) - ave)

                if dimension == 'T':
                    sum /= len(self.T)
                elif dimension == 'S':
                    sum /= len(self.S)
                else:
                    sum /= len(self.G)

                SUM += sum

        if dimension == 'T':
            SUM /= len(self.S) * len(self.G)
        elif dimension == 'S':
            SUM /= len(self.T) * len(self.G)
        else:
            SUM /= len(self.T) * len(self.S)

        return SUM


class Cubes:
    def __init__(self):
        self._cubeVec = list()

    def clear(self):
        for cub in self._cubeVec:
            cub.clear()
        self._cubeVec.clear()

    def contain(self, c: Cube):
        for cub in self._cubeVec:
            if cub.contain(c):
                return True
        return False

    def del_becontained(self, c: Cube):
        to_remove = list()
        for cub in self._cubeVec:
            if c.contain(cub):
                cub.T.clear()
                cub.S.clear()
                cub.G.clear()
                to_remove.append(cub)

        for c in to_remove:
            self._cubeVec.remove(c)

    def addIfMax(self, c: Cube):
        if not self.contain(c):
            self.del_becontained(c)
            self._cubeVec.append(c)

    def validNum(self):
        num = 0
        for cub in self._cubeVec:
            if not cub.deleted:
                num += 1
        return num

    def show(self, array3d):

        sFinal = ""
        i = 0
        sFinal += "===================================================== Clusters ======================================================\n"

        for it in self._cubeVec:
            if not it.deleted:
                sFinal += "\nCluster " + str(i) + ": \n"
                i += 1
                sFinal += "=====================================================================================================================\n"
                sFinal += it.show(array3d)

        return sFinal

    def delet(self, overlaped: float):

        for c in self._cubeVec:
            inc = 0
            total = 0

            for itt in c.T:
                for its in c.S:
                    for itg in c.G:
                        for c2 in self._cubeVec:
                            if c != c2 and not c2.deleted:
                                if c2.containPoint(itt, its, itg):
                                    inc += 1
                                    break
                        total += 1
            if inc/total >= overlaped:
                c.deleted = True

    def merge(self, overlaped: float):
        cont = True

        noMerge = False

        while cont:
            cube = Cube()
            noMerge = True
            it3 = len(self._cubeVec)
            it1 = 0
            while it1 < it3:
                it2 = it1
                it2 += 1
                while it2 < it3:

                    if not self._cubeVec[it1].deleted and not self._cubeVec[it2].deleted:
                        inc = 0
                        total = 0

                        cube.T = SortedSet(
                            sorted(self._cubeVec[it1].T.union(self._cubeVec[it2].T)))
                        cube.S = SortedSet(
                            sorted(self._cubeVec[it1].S.union(self._cubeVec[it2].S)))
                        cube.G = SortedSet(
                            sorted(self._cubeVec[it1].S.union(self._cubeVec[it2].G)))

                        for itt in cube.T:
                            for its in cube.S:
                                for itg in cube.G:
                                    if self._cubeVec[it2].containPoint(itt, its, itg) or self._cubeVec[it1].containPoint(itt, its, itg):
                                        inc += 1
                                    total += 1
                        if inc / total >= overlaped:
                            self._cubeVec[it1].deleted = True
                            self._cubeVec[it2].deleted = True
                            self._cubeVec.append(cube)
                            noMerge = False
                    it2 += 1
                it1 += 1

            if noMerge:
                cont = False

    def showQuality(self, array3d: Array3D):
        tFluc = 0
        sFluc = 0
        gFluc = 0
        cover = dict()

        sFinal = ""

        sFinal += "\n====================================================== Quality ======================================================\n"
        cNum = 0
        eSum = 0
        for it in self._cubeVec:
            if not it.deleted:

                cNum += 1
                eSum += it.elementNum()

                tFluc += it.variance('T', array3d) * it.elementNum()
                sFluc += it.variance('S', array3d) * it.elementNum()
                gFluc += it.variance('G', array3d) * it.elementNum()

                for itt in it.T:
                    for its in it.S:
                        for itg in it.G:
                            cover[itt * (array3d.S * array3d.G) +
                                  its * array3d.G + itg] = 1

        tFluc /= eSum
        sFluc /= eSum
        gFluc /= eSum
        eCover = len(cover)
        cover.clear()

        sFinal += "Clusters#:\t" + str(cNum) + '\n'

        if cNum != 0:
            sFinal += "Elements#:\t" + str(eSum) + '\n'
            sFinal += "Cover:   \t" + str(eCover) + '\n'
            sFinal += "Overlap%:\t{:2.2f}%\n".format(
                float(eSum-eCover) / eCover * 100.0)
            sFinal += "Variance:\tT:{:2.2f}    S:{:2.2f}    G:{:2.2f}\n".format(
                tFluc, sFluc, gFluc)

        return sFinal


class Cluster:

    def __init__(self, array3d: Array3D, cubes: Cubes):
        self._pArray3d = array3d
        self._width = array3d.S
        self._height = array3d.G
        self._edgesMatrix = [list() for _ in range(
            self._width**2)]  # list of list of sets
        self._pRects = [Rects() for _ in range(array3d.T)]  # list of Rects
        self._pCubes = cubes

    def EDGE(self, i, j):
        return self._edgesMatrix[i*self._width+j]

    def EDGE_push(self, i, j, val: set):
        self._edgesMatrix[i*self._width+j].append(val)

    def EDGE_clear(self, i, j):
        self._edgesMatrix[i*self._width+j].clear()

    def getRects(self, T: int, rects: Rects, winsz: float, support: list, eDelta: list, ivn: list):
        """
        ivn: List of float
        eDelta: List of float
        support: List of float
        """
        r = Rect()
        P = SortedSet()
        width = self._pArray3d.S
        edges = 0

        self.getRanges(T, winsz, support[2], ivn)

        # print(self._edgesMatrix)

        i = 0
        while i < width-1:
            j = i + 1
            while j < width:
                edges += len(self.EDGE(i, j))
                j += 1
            i += 1
        print("Time ", str(self._pArray3d.tName[T]), ": Got ", edges,
              " edges, Average edges: ", edges/(self._width*self._width))

        i = 0
        while i < width:
            P.add(i)
            i += 1

        self.EXPAND(rects, r, P, support, eDelta)

        i = 0
        while i < width - 1:
            j = i + 1
            while j < width:
                self.EDGE_clear(i, j)
                j += 1
            i += 1

    def compare_genes(self, g_array, g_rect, ss_rect, winsz, T: int):
        cont = 0
        for s_r in ss_rect:
            val_rect = self._pArray3d.dat(T, s_r, g_rect)
            val_arr = self._pArray3d.dat(T, s_r, g_array)
            if val_rect - winsz <= val_arr <= val_rect + winsz:
                cont += 1
        return cont

    def verify_all_missings(self, rect, T):
        for s in rect.S:
            for g in rect.G:
                if self._pArray3d.dat(T, s, g) != 0:
                    return False
        return True

    def count_mv(self, g_array, ss_rect, T):
        cont = 0
        for s in ss_rect:
            if fZero(self._pArray3d.dat(T, s, g_array)):
                cont += 1
        return cont

    def handle_missings(self, rects: Rects, minG: int, winsz: float, T: int, mthr: float):

        for r in rects.rectVec:
            i = 0
            while i < self._pArray3d.G:
                if i not in r.G and (self.count_mv(i, r.S, T) / len(r.S)) <= mthr:
                    for p_r in r.G:
                        c_v = self.compare_genes(i, p_r, r.S, winsz, T)
                        if c_v == (len(r.S) - self.count_mv(i, r.S, T)):
                            r.G.add(i)
                            break
                i += 1
            if len(r.G) < minG or self.verify_all_missings(r, T):
                r.deleted = True

    def getCubes(self, winsz: float, support: list, eDelta: list, ivn: list, mthr: float):
        c = Cube()
        P = SortedSet()
        width = self._pArray3d.T

        i = 0
        while i < width:
            self.getRects(i, self._pRects[i], winsz, support, eDelta, ivn)

            # CODE ABOUT MISSING DATA ## HERE ## FOR EACH BIC IN TRIC

            self.handle_missings(self._pRects[i], support[0], winsz, i, mthr)

            # print(self._pRects[i].show(self._pArray3d, i))
            print("Biclusters got at this time: ", self._pRects[i].validNum())
            i += 1

        # Ends Biclustering
        i = 0
        while i < width:
            P.add(i)
            i += 1

        self.EXPAND_T(self._pCubes, c, P, support, eDelta)

    def EXPAND_T(self, cubes: Cubes, c: Cube, P: set, support: list, eDelta: list):
        """
        eDelta: List of float
        support: List of float
        """
        minT = support[0]
        minS = support[1]
        minG = support[2]
        c1 = Cube()

        if len(c.T) >= minT and len(c.S) >= minS and len(c.G) >= minG:
            if (eDelta[0] < 0.0 or c.variance('T', self._pArray3d) <= eDelta[0]) and \
                (eDelta[1] < 0.0 or c.variance('S', self._pArray3d) <= eDelta[1]) and \
                    (eDelta[2] < 0.0 or c.variance('G', self._pArray3d) <= eDelta[2]):
                cubes.addIfMax(c)

        while len(P) > 0:
            x = P[0]
            c1.copy(c)
            P.remove(x)
            if len(c1.T) == 0 or c1.T[-1] + 1 == x:  # RESTRIÇÃO TEMPORAL
                c1.T.add(x)
                if len(c.T) == 0:
                    print("Processing Time: ", x)

                for it in self._pRects[x].rectVec:
                    if not it.deleted:
                        if len(c.T) == 0:  # is empty
                            c1.S = copy.deepcopy(it.S)
                            c1.G = copy.deepcopy(it.G)
                        else:
                            c1.S = c.S.intersection(it.S)
                            c1.G = c.G.intersection(it.G)

                        if len(c1.S) >= minS and len(c1.G) >= minG:
                            ctemp = copy.deepcopy(c)
                            self.EXPAND_T(cubes, copy.deepcopy(
                                c1), copy.deepcopy(P), support, eDelta)
                            c.S = ctemp.S
                            c.T = ctemp.T
                            c.G = ctemp.G

    def EXPAND(self, rects: Rects, r: Rect, P, support: list, eDelta: list):
        minS = support[1]
        minG = 1
        r1 = Rect()

        if len(r.S) >= minS and len(r.G) >= minG:
            rects.addIfMax(r)

        while len(P) > 0:  # not empty
            x = P[0]
            r1.copy(r)
            r1.S.add(x)
            P.remove(x)

            if len(r.S) == 0:  # empty
                if not BRIEF_OUT:
                    print("Processing Sample: ", x)
                self.EXPAND(rects, r1, copy.deepcopy(P), support, eDelta)
            else:
                y = r.S[-1]
                for it in self.EDGE(y, x):
                    if len(r.S) == 1:
                        r1.G = SortedSet(sorted(copy.deepcopy(it)))
                    else:
                        r1.G = SortedSet(
                            sorted(r.G.intersection(copy.deepcopy(it))))
                    if len(r1.G) >= minG:
                        rtemp = copy.deepcopy(r)

                        self.EXPAND(rects, copy.deepcopy(r1),
                                    copy.deepcopy(P), support, eDelta)
                        r.G = rtemp.G
                        r.S = rtemp.S

    def getRanges(self, T: int, winsz: float, gSup: int, ivn: list):
        """
        ivn: List of float
        """
        ratioVec = list()  # List of Ratio
        win_size = 1.0 + winsz
        sSet = SortedSet()
        sSet1 = SortedSet()

        i = 0
        while i < self._width - 1:
            j = i + 1
            while j < self._width:
                ratioVec.clear()
                k = 0
                while k < self._height:
                    flag = True
                    # print('i', self._pArray3d.dat(T, i, k))
                    # print('j', self._pArray3d.dat(T, j, k))
                    if not (fZero(self._pArray3d.dat(T, i, k)) or fZero(self._pArray3d.dat(T, j, k))):
                        for iit in ivn:
                            if (fZero(iit - self._pArray3d.dat(T, i, k)) or fZero(iit - self._pArray3d.dat(T, k, k))):
                                flag = False
                                break
                        if flag:
                            ratio = Ratio()

                            ratio.val = round(self._pArray3d.dat(
                                T, j, k) / self._pArray3d.dat(T, i, k), 6)
                            if ratio.val < 0:
                                if self._pArray3d.dat(T, j, k) > 0:
                                    ratio.f1Neg = True
                                else:
                                    ratio.f1Neg = False
                            ratio.trans = k
                            ratioVec.append(ratio)

                    k += 1
                ratioVec.sort()
                vv = ratioVec[0].val
                if vv < 0:
                    for it in ratioVec:
                        it.val -= vv

                while(True):
                    max_sup = sup = 0
                    it1 = it2 = 0

                    while it2 < len(ratioVec):

                        if ratioVec[it1].trans == -1:
                            while it1 < len(ratioVec) and ratioVec[it1].trans == -1:
                                it1 += 1
                            it2 = it1
                        while it2 < len(ratioVec) and ratioVec[it2].trans != -1 and ratioVec[it2].val <= ratioVec[it1].val * win_size:
                            it2 += 1
                            sup += 1

                        if sup > max_sup:
                            start_it = it1
                            end_it = it2
                            max_sup = sup

                        if it2 < len(ratioVec) and ratioVec[it2].trans == -1:
                            it1 = it2
                            sup = 0
                        else:
                            it1 += 1
                            sup -= 1

                    if max_sup >= gSup:
                        it1 = it3 = start_it
                        it2 = it4 = end_it

                        sup = max_sup

                        while sup > gSup:
                            it1 += 1
                            sup -= 1

                        while it2 < len(ratioVec) and ratioVec[it2].trans != -1 and ratioVec[it2].val <= ratioVec[it1].val * win_size:
                            it1 += 1
                            it2 += 1

                        sup = max_sup

                        while sup > gSup:
                            it4 -= 1
                            sup -= 1

                        if it3 != 0:
                            it4 -= 1
                            it3 -= 1
                            while ratioVec[it3].trans != -1 and it3 != 0 and ratioVec[it4].val <= ratioVec[it3].val * win_size:
                                it3 -= 1
                                it4 -= 1
                            if ratioVec[it3].trans == -1 or ratioVec[it4].val > ratioVec[it3].val * win_size:
                                it3 += 1
                        start_it = it3
                        end_it = it2

                        it = start_it
                        while it < end_it:
                            if ratioVec[it].val > 0 or ratioVec[it].f1Neg:
                                sSet.add(ratioVec[it].trans)
                            else:
                                sSet1.add(ratioVec[it].trans)
                            ratioVec[it].trans = -1
                            it += 1

                        if len(sSet) >= gSup:
                            self.EDGE_push(i, j, copy.deepcopy(sSet))

                        if len(sSet1) >= gSup:
                            self.EDGE_push(i, j, copy.deepcopy(sSet1))

                        sSet.clear()
                        sSet1.clear()

                    if max_sup < gSup:
                        break
                ratioVec.clear()
                j += 1
            i += 1


if __name__ == "__main__":
    import argparse
    from sys import argv
    array3d = None
    cubes = None

    support = [None] * 3  # INT       0:T, 1:S, 2:G
    showClusters = [False] * 3  # BOOL  0:T, 1:S, 2:G
    eDelta = [0] * 3  # FLOAT       0:T, 1:S, 2:G
    invalidNum = list()

    win_size = 0.03
    del_overlap = 1.00
    mer_overlap = 1.00
    mv_threshold = 0.0

    ###########################################################################################################

    parser = argparse.ArgumentParser(description='TCtriCluster')
    parser.add_argument('-f', action='store', dest='filename',
                        type=str, help='/*File name', required=True)

    parser.add_argument('-sT', action='store', dest='minT', type=int,
                        help='/* Minimum size for T dimension', required=True)
    parser.add_argument('-sS', action='store', dest='minS', type=int,
                        help='/* Minimum size for S dimension', required=True)
    parser.add_argument('-sG', action='store', dest='minG', type=int,
                        help='/* Minimum size for G dimension', required=True)

    parser.add_argument('-w', action='store', dest='winsz',
                        type=float, help='Range Window Size: /*0.03 by default*/')
    parser.add_argument('-d', action='store', dest='dele',
                        type=float, help='Deletion Threshold   /*1.00 by default*/')
    parser.add_argument('-m', action='store', dest='merg', type=float,
                        help='Merging  Threshold:    /*1.00 by default*/')
    parser.add_argument('-mv', action='store', dest='mv', type=float,
                        help='Missing Values  Threshold:    /*0.00 by default*/')

    parser.add_argument('-u', action='store', dest='unrel', type=float,
                        help='Unrelated Numbers: /*mciroCluster will not consider the values refered by this option*/')

    parser.add_argument('-et', action='store', dest='et', type=float,
                        help='when trying to get close Time values only')
    parser.add_argument('-es', action='store', dest='es', type=float,
                        help='when trying to get close Sample values only')
    parser.add_argument('-eg', action='store', dest='eg', type=float,
                        help='when trying to get close Gene values only')

    parser.add_argument('-r', action='store', dest='recfile',
                        type=str, help='Output the result to a file')
    parser.add_argument('-b', action='store', dest='brief', type=bool,
                        default=False, help='/*Output the current status in brief*/')
    parser.add_argument('-o', action='store', dest='option', type=int,
                        help='1:Original clusters,  2:Clusters  after deletion, 3:Clusters  after merging*/')

    for i in range(3):
        showClusters[i] = False
        eDelta[i] = -1.0

    results = parser.parse_args()
    filename = results.filename
    support[0] = results.minT
    support[1] = results.minS
    support[2] = results.minG

    if results.winsz is not None:
        win_size = results.winsz
    if results.dele is not None:
        del_overlap = results.dele
    if results.merg is not None:
        mer_overlap = results.merg
    if results.mv is not None:
        mv_threshold = results.mv

    if results.unrel is not None:
        invalidNum.append(results.unrel)

    if results.brief:
        BRIEF_OUT = True

    if results.et is not None:
        eDelta[0] = results.et
    if results.es is not None:
        eDelta[1] = results.es
    if results.eg is not None:
        eDelta[2] = results.eg

    if results.option == 1:
        showClusters[0] = True
    elif results.option == 2:
        showClusters[1] = True
    elif results.option == 3:
        showClusters[2] = True

    print(" ".join(argv)+'\n')

    print("Input File Name:\t", filename)
    print("Cluster Minimum Size:\tT:",
          support[0], ", S:", support[1], ", G:", support[2])
    print("Ratio Window Size:\t", win_size)
    print("Deletion Threshold:\t", del_overlap)
    print("Merging  Threshold:\t", mer_overlap)
    print("Missing Values  Threshold:\t", mv_threshold)

    if not fZero(eDelta[0]+1.0):
        print("delta-T:\t\t", eDelta[0])
    if not fZero(eDelta[1]+1.0):
        print("delta-S:\t\t", eDelta[1])
    if not fZero(eDelta[2]+1.0):
        print("delta-G:\t\t", eDelta[2])

    print()
    print()

    ###########################################################################################################

    cubes = Cubes()
    inp = Input(filename)
    array3d = inp.readTabFile()
    cluster = Cluster(array3d, cubes)

    cluster.getCubes(win_size, support, eDelta, invalidNum, mv_threshold)
    print()
    print("Original clusters: ", cubes.validNum())

    if showClusters[0]:
        print(cubes.show(cluster._pArray3d))

    print(cubes.showQuality(array3d))

    if not fZero(del_overlap - 1.0):
        cubes.delet(del_overlap)
        print("Clusters after deletion: ", cubes.validNum())
        if showClusters[1]:
            print(cubes.show(array3d))
        print(cubes.showQuality(array3d))

    if not fZero(mer_overlap - 1.0):
        cubes.merge(mer_overlap)
        print("Clusters after merging: ", cubes.validNum())
        if showClusters[2]:
            print(cubes.show(array3d))
        print(cubes.showQuality(array3d))
