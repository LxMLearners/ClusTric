from tricluster import Tricluster
import numpy as np
import pandas as pd
import math
import scipy.stats as ss
import statistics as stat


def get_triclusters(output_file):
    fp = open(output_file, 'r')
    content = fp.readlines()

    content = list(map(lambda x: x.strip(), content))
    flag = False
    triclusters = list()
    for c in content:
        if c.startswith("====================================================================================================================="):
            flag = True
        elif c.startswith("|T|x|S|x|G|") and flag:
            l = c.split(":")[1].split('x')
            tri = Tricluster(int(l[0]), int(l[1]), int(l[2]))
        elif c.startswith("Time") and flag:
            time = int(c.split(":")[1].strip())
            tri.addTime("T-" + str(time))
        elif c.startswith("S-") and flag:
            smp = list(map(lambda x: int(x.strip()), c.split("S-")[1:]))
            for s in smp:
                tri.addSample("S-" + str(s))
        elif c.startswith("G") and flag:
            pc = c.split('\t')[0]
            tri.addPatient(pc)
            q = c.split("\t")[1].strip().replace(" ", "$")
            vals = q.split("$")
            vals = list(filter(lambda z: len(z) > 0, vals))
            idx = 0
            for v in vals:
                s = tri.getSamples()[idx]
                tri.addValue("T-" + str(time), s, pc, float(v.strip()))
                idx += 1
        elif c.startswith("Cluster") and flag:
            triclusters.append(tri)
    fp.close()
    return triclusters


def compute_representative_patterns(tricluster, categorical_feats, continuos_feats):
    """
    Format: [[], []]
    """

    result = []
    for tp in tricluster.getTimes():
        bic_vals = list()
        bicluster = [tricluster.getFeatValues(
            tp, f) for f in tricluster.getSamples()]

        means = list(map(lambda z: aux_mean(z), bicluster))
        modes = list(map(lambda z: aux_mode(z), bicluster))

        i = 0
        for ft in tricluster.getSamples():
            if ft in categorical_feats:
                bic_vals.append((ft, modes[i]))
            else:
                bic_vals.append((ft, means[i]))
            i += 1
        result.append(bic_vals)
    return(result)


def virtual_3d(bics_list, tricluster, categorical_feats, continuos_feats):
    # bics_list is a list of list representing the biclusters 2D v. pattern
    v_pat_3D = list()
    for i in range(len(bics_list[0])):
        comp_l = list()
        for bic_p in bics_list:
            comp_l.append(bic_p[i])

        if tricluster.getSamples()[i] in categorical_feats:
            p = aux_mode(comp_l)
            if type(p) != np.float64:
                p = list(p)
            v_pat_3D.append(p)
        else:
            if type(comp_l[0]) == tuple:
                c = list(map(lambda a: a[1], comp_l))
                p = aux_mean(c)
                p = [comp_l[0][0], p]
            else:
                p = aux_mean(comp_l)
            v_pat_3D.append(p)
    return v_pat_3D


def representative_pattern_3d(tricluster, categorical_feats, continuos_feats):

    bics_list = compute_representative_patterns(
        tricluster, categorical_feats, continuos_feats)
    return(virtual_3d(bics_list, tricluster, categorical_feats, continuos_feats))


def compute_similarity_matrix(patients_tri, triclusters, categorical_feats, continuos_feats, corr=False, tri=False):
    """

    Parameters
    ----------
    patients_tri: list of Tricluster object with patients information, 
    triclusters: list of Triclusters
    categorical_feats: list of categorical features, 
    continuos_feats: list of continuos features, 
    corr: `True`if uses Pearson correlation as metric. `False` for Euclidean Distance 
    tri: `True` if similarities are computed using triclusters, `False` (default) if using biclusters
    """

    dist_set = dict()
    print(len(patients_tri))

    p_pats_reps = list()
    for bic_p in patients_tri:
        # Compute Patient Patterns
        p_pats = [bic_p.getPatientsVals(t=tp) for tp in bic_p.getTimes()]
        if tri:
            p_pats = p_pats[:-1]
            p_pats = virtual_3d(
                p_pats, bic_p, categorical_feats, continuos_feats)
        p_pats_reps.append(p_pats)

    tric_names = list()
    t_pats_reps = list()
    i = 0
    for tric in triclusters:  # Computing Triclusters Virtual Pattern
        if tri:
            t_pats = representative_pattern_3d(
                tric, categorical_feats, continuos_feats)
            if t_pats not in t_pats_reps:
                t_pats_reps.append(t_pats)
                tric_names.append("Tric_" + str(i))
        else:
            t_pats = compute_representative_patterns(
                tric, categorical_feats, continuos_feats)  # [[Pattern T0],[Pattern T1],...]
            if t_pats not in t_pats_reps:
                t_pats_reps.append(t_pats)
                tric_names.append("Tric_" + str(i))
        i += 1

    pat = 0
    for p_pats in p_pats_reps:
        nTr = 0
        for t in t_pats_reps:
            if tri:
                # features indices
                idx = list(map(lambda x: int(x[0].split('-')[-1]), t))
                # pattern tricluster
                a = np.array(tuple(list(map(lambda y: float(y[1]), t))))
                b = np.array(tuple(p_pats))[[idx]]  # pattern patient

                if not corr:
                    dist_set[(pat, tric_names[nTr])] = np.linalg.norm(a-b)
                else:
                    if len(a) == 2:
                        dist_set[(pat, tric_names[nTr])] = dot_product(
                            list(a), list(b[0]))
                    else:
                        dist_set[(pat, tric_names[nTr])] = pearson_correlation(
                            list(a), list(b[0]))
            else:
                idx = list(map(lambda x: int(x[0].split('-')[-1]), t[0]))
                for i in range(len(t)):

                    # pattern tricluster
                    a = np.array(tuple(list(map(lambda y: y[1], t[i]))))
                    b = np.array(tuple(p_pats[i]))[[idx]]  # pattern patient

                    if not corr:
                        dist_set[(pat, tric_names[nTr], i)
                                 ] = np.linalg.norm(a-b)
                    else:
                        if len(a) == 2:
                            dist_set[(pat, tric_names[nTr], i)
                                     ] = dot_product(list(a), list(b[0]))
                        else:
                            dist_set[(pat, tric_names[nTr], i)
                                     ] = pearson_correlation(list(a), list(b[0]))
            nTr += 1
        pat += 1

    f_matrix = list()
    for p in range(len(patients_tri)):
        vals = filter(lambda x: x[0] == p, dist_set.keys())
        line = [dist_set[v] for v in vals]
        f_matrix.append(line)

    vals = filter(lambda x: x[0] == 0, dist_set.keys())
    if tri:
        cols = list(map(lambda z: str(z[1]), vals))
    else:
        cols = list(map(lambda z: str(z[1]) + "_" + str(z[2]), vals))

    return f_matrix, cols


def write_matrix(final_matrix, file_name, classes, cols, prog_r=None):
    fin = pd.DataFrame(data=final_matrix, columns=cols)
    fin["Evolution"] = classes
    if prog_r is not None:
        fin["ProgRate"] = prog_r
    fin.to_csv(file_name, index=False)

######## AUXILIARES ##########


def aux_mean(z):
    z = list(filter(lambda x: x != 0, z))
    if len(z) == 0:
        return 0
    else:
        return stat.mean(z)


def aux_mode(z):
    if len(z) == 0:
        return 0
    else:
        return ss.mode(z, keepdims = True)[0][0]


def pearson_correlation(numbers_x, numbers_y):

    for a, b in list(zip(numbers_x, numbers_y)):
        if a == math.nan or b == math.nan:
            return math.nan

    mean_x = sum(numbers_x) / len(numbers_x)
    mean_y = sum(numbers_y) / len(numbers_y)

    subtracted_mean_x = [i - mean_x for i in numbers_x]
    subtracted_mean_y = [i - mean_y for i in numbers_y]

    subtracted_mean_square_x = [(i - mean_x)**2 for i in numbers_x]
    subtracted_mean_square_y = [(i - mean_y)**2 for i in numbers_y]

    x_times_y = [a * b for a,
                 b in list(zip(subtracted_mean_x, subtracted_mean_y))]

    dem = (math.sqrt(sum(subtracted_mean_square_x))
           * math.sqrt(sum(subtracted_mean_square_y)))

    if dem == 0:
        return 0

    result = sum(x_times_y) / dem

    return result


def dot_product(lst1, lst2):

    a1, a2 = lst1
    b1, b2 = lst2

    if a1 == math.nan or a2 == math.nan or b1 == math.nan or b2 == math.nan:
        return math.nan

    inner_product = a1*b1 + a2*b2

    len1 = math.hypot(a1, a2)
    len2 = math.hypot(b1, b2)

    dem = len1*len2
    if dem == 0:
        return 0

    result = inner_product/dem

    return result
