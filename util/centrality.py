import time
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import copy


def unit_vector(vector):
    """Returns the unit vector of the vector."""
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2):
    """Returns the angle in radians between vectors 'v1' and 'v2'::

    >>> angle_between((1, 0, 0), (0, 1, 0))
    1.5707963267948966
    >>> angle_between((1, 0, 0), (1, 0, 0))
    0.0
    >>> angle_between((1, 0, 0), (-1, 0, 0))
    3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


# MultiCens: local centrality
def local_centrality(A_tilde_full, num_layers, p):
    print("\n[Computing Local Centrality]\n")
    print("Im in")
    print(A_tilde_full)
    print("im something")
    print(num_layers)
    print("By")
    print(p)
    start_time = time.time()
    A_tilde_full = A_tilde_full / np.sum(A_tilde_full, axis=0)
    # A_tilde_full = A_tilde_full/LA.norm(A_tilde_full)
    n = int(np.shape(A_tilde_full)[0] / num_layers)
    N = int(np.shape(A_tilde_full)[0])
    A_tilde = np.zeros_like(A_tilde_full, dtype=np.float32)
    print("middle")
    for i in range(num_layers):
        A_tilde[(i * n) : ((i + 1) * n), (i * n) : ((i + 1) * n)] = A_tilde_full[
            (i * n) : ((i + 1) * n), (i * n) : ((i + 1) * n)
        ]

    ones_t = np.ones((N,)) / N
    l = np.copy(ones_t)
    l_new = np.copy(ones_t)

    

    count = 0
    current_angle = np.zeros(
        3,
    )
    # print("Local Centrality computation starts")
    print("somewhere")
    while count < 200:
        # print(count)
        count = count + 1
        l_new = (p * ((A_tilde).dot(l))) + (1 - p) * ones_t
        current_angle[0] = current_angle[1]
        current_angle[1] = current_angle[2]
        current_angle[2] = angle_between(l, l_new)
        # print(current_angle)
        if ((current_angle[1] == current_angle[2])) or (current_angle[2] == 0):
            break
        l = copy.deepcopy(l_new)
        # l = copy.deepcopy(l_new)/LA.norm(l_new)
    # print("l values")
    # print(l)
    # new_l = l/l.sum()# + l_tt/l_tt.sum()
    new_l = copy.deepcopy(l)
    # ---new_l is empty
    for i in range(num_layers):
        l[(i * n) : ((i + 1) * n)] = (
            l[(i * n) : ((i + 1) * n)] / l[(i * n) : ((i + 1) * n)].sum()
        )
    
    new_l = new_l / np.max(new_l)
    print("end")
    # Time calculation
    end_time = time.time()
    execution_time = end_time - start_time
    execution_time = round(execution_time/60, 2)
    time_taken = f'{execution_time} min(s)'
    print("\nTime taken to compute the local centrality scores:", time_taken, "\n")
    return new_l


# MultiCens: global centrality
def global_centrality(A_tilde_full, num_layers, p):
    print("\n[Computing Global Centrality]\n")
    start_time = time.time()
    l = local_centrality(A_tilde_full, num_layers, p)  # last checkpoint
    # A_tilde = copy.deepcopy(A_tilde_full)
    A_tilde = A_tilde_full / np.sum(A_tilde_full, axis=0)
    N = int(np.shape(A_tilde)[0])
    n = int(N / num_layers)
    # print("A_tilde sums")
    # print(np.sum(A_tilde[0,:]))
    # print(np.sum(A_tilde[:,0]))

    A = np.zeros_like(A_tilde_full, dtype=np.float32)
    C = np.zeros_like(A_tilde_full, dtype=np.float32)
    A[:n, :n] = copy.deepcopy(A_tilde[:n, :n])
    A[n : (2 * n), n : (2 * n)] = copy.deepcopy(A_tilde[n : (2 * n), n : (2 * n)])
    A[(2 * n) : (3 * n), (2 * n) : (3 * n)] = copy.deepcopy(
        A_tilde[(2 * n) : (3 * n), (2 * n) : (3 * n)]
    )
    A[(3 * n) : (4 * n), (3 * n) : (4 * n)] = copy.deepcopy(
        A_tilde[(3 * n) : (4 * n), (3 * n) : (4 * n)]
    )
    C = A_tilde - A

    ones_t = np.ones((N,)) / N
    g = copy.deepcopy(ones_t)
    g_new = copy.deepcopy(ones_t)

    counter = 0
    current_angle = np.zeros(
        3,
    )

    while counter < 150:
        g_new = (p * ((A + C).dot(g) + (C.dot(l)))) + ((1 - p) * ones_t)

        current_angle[0] = current_angle[1]
        current_angle[1] = current_angle[2]
        current_angle[2] = angle_between(g, g_new)
        # print(current_angle)
        if ((current_angle[1] == current_angle[2])) or (current_angle[2] == 0):
            break
        g = copy.deepcopy(g_new)
        # g = copy.deepcopy(g_new)/LA.norm(g_new)
        counter += 1
        # print(counter)
        # print(np.sum(g))

    new_g = copy.deepcopy(g)
    g_fresh = copy.deepcopy(g)
    new_g = new_g / new_g.sum()
    g_fresh = g_fresh / g_fresh.sum()
    
    new_g = new_g / np.max(new_g)
    
    # Time calculation
    end_time = time.time()
    execution_time = end_time - start_time
    execution_time = round(execution_time/60, 2)
    time_taken = f'{execution_time} min(s)'
    print("\nTime taken to compute the global centrality scores:", time_taken, "\n")
    
    return new_g


def right_new_local_centrality_st(
    A_tilde_full, num_layers, target_tissue, target_gene_indices, start, end, p
):
    # print("target tissue id")
    # print(target_tissue)
    # print("right_new_local_centrality")
    A_tilde_full = A_tilde_full / np.sum(A_tilde_full, axis=0)
    # A_tilde_full = A_tilde_full/LA.norm(A_tilde_full)
    num_target_genes = int(np.shape(target_gene_indices)[0])
    n = int(np.shape(A_tilde_full)[0] / num_layers)
    N = int(np.shape(A_tilde_full)[0])
    A_tilde = np.zeros_like(A_tilde_full, dtype=np.float32)

    # if target_tissue == 1:
    #     A_tilde[n:,n:] = A_tilde_full[n:,n:]
    # elif target_tissue == 0:
    #     A_tilde[:n,:n] = A_tilde_full[:n,:n]
    # else:
    #     print("wrong target tissue")

    A_tilde[start:end, start:end] = A_tilde_full[
        start:end, start:end
    ]  # replacement for above commented code

    np.fill_diagonal(A_tilde, 0.0)
    # A_tilde = A_tilde/LA.norm(A_tilde)

    ones_t = np.zeros((N,))
    # ones_t[
    #     (np.asarray(target_gene_indices, dtype=np.int32) + int(n * target_tissue))
    # ] = 1 / len(target_gene_indices)

    ones_t[np.asarray(target_gene_indices, dtype=np.int32)] = 1 / len(
        target_gene_indices
    )  # replacement for above commented code

    l = np.copy(ones_t)
    l_new = np.copy(ones_t)

    count = 0
    current_angle = np.zeros(
        3,
    )
    print("[Local Centrality computation starts]")
    while count < 200:
        # print(count)
        count = count + 1
        l_new = (p * ((A_tilde).dot(l))) + (1 - p) * ones_t
        current_angle[0] = current_angle[1]
        current_angle[1] = current_angle[2]
        current_angle[2] = angle_between(l, l_new)
        # print(current_angle)
        if (
            (current_angle[0] == current_angle[1])
            and (current_angle[1] == current_angle[2])
        ) or (current_angle[2] == 0):
            break
        l = copy.deepcopy(l_new)
        # l = copy.deepcopy(l_new)/LA.norm(l_new)
    # print("l values")
    # print(l)
    # new_l = l/l.sum()# + l_tt/l_tt.sum()
    new_l = copy.deepcopy(l)
    # if target_tissue == 0:
    #     print("Found local centrality for target set centrality")
    #     new_l[:n] = l[:n] / l[:n].sum()
    # elif target_tissue == 1:
    #     new_l[n:] = l[n:] / l[n:].sum()
    #     print("Found local centrality for source set centrality")
    # else:
    #     print("invalid target tissue")

    new_l[start:end] = (
        l[start:end] / l[start:end].sum()
    )  # replacement for above commented code

    # print("Target genes")
    # print(target_gene_indices)
    return new_l


def right_target_global_centrality_t(
    A_tilde_full, num_layers, target_tissue, target_gene_indices, start, end, p
):
    print("\n[Computing Query-set Centrality]\n")
    start_time = time.time()
    l = right_new_local_centrality_st(
        A_tilde_full, num_layers, target_tissue, target_gene_indices, start, end, p
    )  # last checkpoint
    # A_tilde = copy.deepcopy(A_tilde_full)
    A_tilde = A_tilde_full / np.sum(A_tilde_full, axis=0)
    N = int(np.shape(A_tilde)[0]) # rows in correlation matrix
    n = int(N / num_layers) # number of genes selected from each tissue
    # print("A_tilde sums")
    # print(np.sum(A_tilde[0, :]))
    # print(np.sum(A_tilde[:, 0]))

    A = np.zeros_like(A_tilde_full, dtype=np.float32)
    C = np.zeros_like(A_tilde_full, dtype=np.float32)
    # A[:n, :n] = copy.deepcopy(A_tilde[:n, :n])
    # A[n:, n:] = copy.deepcopy(A_tilde[n:, n:])
    # C[:n, n:] = copy.deepcopy(A_tilde[:n, n:])
    # C[n:, :n] = copy.deepcopy(A_tilde[n:, :n])
    
    # replacement for above commented code
    # N is the number of sub-matrices
    # n size of each sub-matrix
    for i in range(num_layers):
        for j in range(num_layers):
            if i == j:
                A[i*n:(i+1)*n, j*n:(j+1)*n] = copy.deepcopy(A_tilde[i*n:(i+1)*n, j*n:(j+1)*n])
            else:
                C[i*n:(i+1)*n, j*n:(j+1)*n] = copy.deepcopy(A_tilde[i*n:(i+1)*n, j*n:(j+1)*n])

    print(f"[Finding target global centrality for layer {str(target_tissue+1)}]")

    ones_t = np.zeros((N,))
    # ones_t[
    #     (np.asarray(target_gene_indices, dtype=np.int32) + int(n * target_tissue))
    # ] = 1 / len(target_gene_indices)

    ones_t[np.asarray(target_gene_indices, dtype=np.int32)] = 1 / len(
        target_gene_indices
    )  # replacement for above commented code
    g = copy.deepcopy(ones_t)
    g_new = copy.deepcopy(ones_t)

    counter = 0
    current_angle = np.zeros(
        3,
    )
    while counter < 200:
        g_new = (p * ((A + C).dot(g) + (C.dot(l)))) + ((1 - p) * ones_t)

        current_angle[0] = current_angle[1]
        current_angle[1] = current_angle[2]
        current_angle[2] = angle_between(g, g_new)
        # print(current_angle)
        if (
            (current_angle[0] == current_angle[1])
            and (current_angle[1] == current_angle[2])
        ) or (current_angle[2] == 0):
            break
        g = copy.deepcopy(g_new)
        # g = copy.deepcopy(g_new)/LA.norm(g_new)
        counter += 1
        # print(counter)
        # print(np.sum(g))
    new_g = copy.deepcopy(g)
    # if target_tissue == 1:
    #     print("Found target set centrality")
    #     new_g[:n] = g[:n] / g[:n].max()
    # elif target_tissue == 0:
    #     new_g[n:] = g[n:] / g[n:].max()
    #     print("Found source set centrality")
    # else:
    #     print("invalid target tissue")

    new_g[start:end] = (
        g[start:end] / g[start:end].max()
    )  # replacement for above commented code
    
     # Time calculation
    end_time = time.time()
    execution_time = end_time - start_time
    execution_time = round(execution_time/60, 2)
    time_taken = f'{execution_time} min(s)'
    print("\nTime taken to compute the query-set centrality scores:", time_taken, "\n")
    
    return l, new_g


def plot_k_curve(genes, cen_vector, ground_truth_genes, filtered, n):
    secreted_proteins = list(
        pd.read_csv("../data/002790_proteins.csv", header=None).values
    )
    secreted_proteins = [s[0].upper() for s in secreted_proteins]
    order = cen_vector[:n].argsort()
    ranks = order.argsort()
    d = {"gene_name": genes, "centrality": cen_vector[:n], "rank": n - ranks}
    results = pd.DataFrame(data=d)
    results = results.sort_values(by=["centrality"], ascending=False)
    ranked_genes = results["gene_name"].tolist()

    filtered_results = results[results.gene_name.isin(secreted_proteins)]
    filtered_ranked_genes = filtered_results["gene_name"].tolist()

    if filtered:
        n_filtered = filtered_results.shape[0]
        k_range = np.arange(n_filtered)
        hits = np.zeros_like(k_range)
        current_hit_count = 0
        for i in k_range:
            if filtered_ranked_genes[i] in ground_truth_genes:
                current_hit_count = current_hit_count + 1
            hits[i] = current_hit_count
        random_curve = np.cumsum(np.full((n_filtered,), hits[-1] / n_filtered))
        plot_variables = {"recall_at_k": hits, "random_curve": random_curve}
        plot_df = pd.DataFrame(data=plot_variables)
        sns.lineplot(data=plot_df)

    else:
        k_range = np.arange(n)
        hits = np.zeros_like(k_range)
        current_hit_count = 0
        for i in k_range:
            if ranked_genes[i] in ground_truth_genes:
                current_hit_count = current_hit_count + 1
            hits[i] = current_hit_count
        random_curve = np.cumsum(np.full((n,), hits[-1] / n))
        plot_variables = {"query-set centrality": hits, "random_curve": random_curve}
        plot_df = pd.DataFrame(data=plot_variables)
        sns.lineplot(data=plot_df)
        # plt.xlabel("Top k predictions", labelpad=1)
        # plt.ylabel("Recall at k", labelpad=1)
        plt.xticks(np.arange(0, 15001, 3500))
        plt.savefig(
            "./insulin_responding_results_corr_SNAP.svg",
            dpi=300,
            fontsize=17,
            bbox_inches="tight",
        )
        # plt.title("My Daily Step Count Tracked by Fitbit", y=1.02, fontsize=22);
        max_area = (current_hit_count * (current_hit_count + 1)) / 2 + (
            (n - current_hit_count) * current_hit_count
        )
        method_area = np.sum(hits)
        print("Area under curve is: ", (method_area / max_area))

    # get lncRNAs ranking
    lncRNAs = list(np.load("../data/paper/lncRNAs.npy"))
    lncRNA_results = results[results.gene_name.isin(lncRNAs)]
    lncRNA_ranked_genes = lncRNA_results["gene_name"].tolist()

    print(hits[-1])
    return (plot_df, results, filtered_results, lncRNA_results)
