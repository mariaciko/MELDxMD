#!/usr/bin/env python
# encoding: utf-8

import argparse

import numpy as np
import progressbar
from matplotlib import pyplot

from meld import vault


def main():
    store = vault.DataStore.load_data_store()
    store.initialize(mode="r")

    parser = argparse.ArgumentParser(
        description="Analyze the results of replica exchange."
    )

    parser.add_argument(
        "--start", type=int, default=None, help="first frame to extract (default: first"
    )
    parser.add_argument(
        "--end", type=int, default=None, help="last frame to extract (default: last)"
    )
    subparsers = parser.add_subparsers(dest="command")

    viz_alpha = subparsers.add_parser(
        "visualize_alpha", help="visualize the value of alpha"
    )
    viz_alpha.set_defaults(func=visualize_alpha)

    ext_alpha = subparsers.add_parser("extract_alpha", help="extract alpha values")
    ext_alpha.add_argument("outfile", help="output filename")
    ext_alpha.set_defaults(func=extract_alpha)

    viz_trace = subparsers.add_parser(
        "visualize_trace", help="visualize trace of replica through ladder"
    )
    viz_trace.add_argument(
        "--replicas", nargs="+", type=int, default=None, help="replicas to visualize"
    )
    viz_trace.set_defaults(func=visualize_trace)

    ext_trace = subparsers.add_parser("extract_trace", help="extract trace")
    ext_trace.add_argument("outfile", help="output filename")
    ext_trace.set_defaults(func=extract_trace)

    viz_fup = subparsers.add_parser("visualize_fup", help="visualize the value of f_up")
    viz_fup.set_defaults(func=visualize_fup)

    ext_fup = subparsers.add_parser("extract_fup", help="extract fup")
    ext_fup.add_argument("outfile", help="output filename")
    ext_fup.set_defaults(func=extract_fup)

    viz_accept = subparsers.add_parser(
        "visualize_accept", help="visualize acceptance probabilities"
    )
    viz_accept.set_defaults(func=visualize_accept)

    viz_alpha_vs_rep = subparsers.add_parser(
        "visualize_alpha_vs_rep", help="visualize alpha versus replica index"
    )
    viz_alpha_vs_rep.add_argument("--frame", type=int, help="frame to vizualize")
    viz_alpha_vs_rep.set_defaults(func=alpha_vs_rep)

    viz_accept_vs_rep = subparsers.add_parser(
        "visualize_accept_vs_rep", help="visualize acceptance rate versus replica index"
    )
    viz_accept_vs_rep.add_argument("--frame", type=int, help="frame to vizualize")
    viz_accept_vs_rep.set_defaults(func=accept_vs_rep)

    args = parser.parse_args()
    args.func(store, args)


def alpha_vs_rep(store, args):
    frame = args.frame
    alphas = store.load_alphas(frame)
    pyplot.plot(alphas, "ko-")
    pyplot.show()


def accept_vs_rep(store, args):
    frame = args.frame
    alphas = store.load_acceptance_probabilities(frame)
    pyplot.plot(alphas, "ko-")
    pyplot.show()


def visualize_alpha(store, args):
    alphas = get_alphas(store, args.start, args.end)
    n_reps = alphas.shape[0]
    for i in range(n_reps):
        pyplot.plot(alphas[i, :], "k-")
    pyplot.show()


def extract_alpha(store, args):
    alphas = get_alphas(store, args.start, args.end)
    np.savetxt(args.outfile, alphas)


def visualize_trace(store, args):
    perm_vecs = get_permutation_vectors(store, args.start, args.end)
    traces = deshuffle_traces(perm_vecs)

    n_replicas = traces.shape[1]
    n_steps = traces.shape[0]

    current_index = np.array(list(range(n_replicas)))

    results = []
    for step in range(n_steps):
        new_value = np.zeros_like(current_index)
        new_value[traces[step, :]] = current_index
        results.append(new_value)
    results = np.array(results)

    if args.replicas is None:
        reps = range(n_replicas)
    else:
        reps = args.replicas

    for index in reps:
        pyplot.plot(results[:, index])

    pyplot.show()


def extract_trace(store, args):
    perm_vecs = get_permutation_vectors(store, args.start, args.end)
    traces = deshuffle_traces(perm_vecs)
    np.savetxt(args.outfile, traces)


def roundtrips(store, args):
    perm_vecs = get_permutation_vectors(store, args.start, args.end)
    traces = deshuffle_traces(perm_vecs)
    n_reps = traces.shape[1]


def visualize_fup(store, args):
    perm_vecs = get_permutation_vectors(store, args.start, args.end)
    traces = deshuffle_traces(perm_vecs)
    n_reps = traces.shape[1]
    f_up, f_down = compute_fup(traces)
    pyplot.plot(f_up, "k-")
    pyplot.plot([0, n_reps - 1], [1, 0], "k--")
    pyplot.plot(f_down, "r-")
    pyplot.plot([0, n_reps - 1], [0, 1], "r--")
    pyplot.xlim(0, n_reps - 1)
    pyplot.show()


def visualize_accept(store, args):
    if args.start:
        start = args.start
    else:
        start = 1

    if args.end:
        end = args.end
    else:
        end = store.max_safe_frame

    n_steps = end - start
    n_pairs = store.n_replicas - 1

    bar = get_progress_bar("Loading data", n_steps).start()
    probs = np.zeros((n_pairs, n_steps))
    for index, frame in enumerate(range(start, end)):
        bar.update(index)
        probs[:, index] = store.load_acceptance_probabilities(frame)
    bar.finish()

    for i in range(n_pairs):
        pyplot.plot(probs[i, :])
    pyplot.show()


def extract_fup(store, args):
    perm_vecs = get_permutation_vectors(store, args.start, args.end)
    traces = deshuffle_traces(perm_vecs)
    f_up, f_down = compute_fup(traces)
    np.savetxt(args.outfile, f_up)


def get_alphas(store, start, end):
    if start is None:
        start = 1
    if end is None:
        end = store.max_safe_frame

    bar = get_progress_bar("Loading data", end - start).start()

    alphas = np.zeros((store.n_replicas, end - start))
    for index, frame in enumerate(range(start, end)):
        bar.update(index)
        alphas[:, index] = store.load_alphas(frame)
    bar.finish()

    return alphas


def get_permutation_vectors(store, start, end):
    if start is None:
        start = 1
    if end is None:
        end = store.max_safe_frame - 1

    perm_vecs = np.zeros((store.n_replicas, end - start), dtype=int)
    bar = get_progress_bar("Loading data", end - start).start()

    for index, frame in enumerate(range(start, end)):
        bar.update(index)
        perm_vecs[:, index] = store.load_permutation_vector(frame)
    bar.finish()

    return perm_vecs


def deshuffle_traces(perm_vecs):
    n_replicas = perm_vecs.shape[0]
    n_steps = perm_vecs.shape[1]

    results = []
    current_indices = np.array(list(range(n_replicas)))

    bar = get_progress_bar("Deshuffling", n_steps).start()
    for i in range(n_steps):
        bar.update(i)
        current_indices = current_indices[perm_vecs[:, i]]
        results.append(current_indices)
    bar.finish()
    return np.array(results)


def calc_roundtrips(traces):
    # these are indexed by LADDER STEP
    up_count = np.zeros_like(traces[0, :])
    down_count = np.zeros_like(traces[0, :])
    # these are indexed by REPLICA
    going_up = np.zeros_like(traces[0, :])
    going_down = np.zeros_like(traces[0, :])

    n_steps = traces.shape[0]
    n_reps = traces.shape[1]

    bar = get_progress_bar("Computing f_up", n_steps).start()
    for step in range(n_steps):
        bar.update(step)
        # this gives the REPLICA for each LADDER STEP
        this_trace = traces[step, :]

        # these are the top and bottom REPLICAS
        top = this_trace[n_reps - 1]
        bottom = this_trace[0]

        # update if each replica is going up or down
        going_up[top] = 0
        going_down[top] = 1
        going_up[bottom] = 1
        going_down[bottom] = 0

        up_count += going_up[this_trace]
        down_count += going_down[this_trace]
    bar.finish()

    f_up = up_count / (up_count + down_count + 0.1)
    f_down = down_count / (up_count + down_count + 0.1)
    return f_up, f_down


def compute_fup(traces):
    # these are indexed by LADDER STEP
    up_count = np.zeros_like(traces[0, :])
    down_count = np.zeros_like(traces[0, :])
    # these are indexed by REPLICA
    going_up = np.zeros_like(traces[0, :])
    going_down = np.zeros_like(traces[0, :])

    n_steps = traces.shape[0]
    n_reps = traces.shape[1]

    bar = get_progress_bar("Computing f_up", n_steps).start()
    for step in range(n_steps):
        bar.update(step)
        # this gives the REPLICA for each LADDER STEP
        this_trace = traces[step, :]

        # these are the top and bottom REPLICAS
        top = this_trace[n_reps - 1]
        bottom = this_trace[0]

        # update if each replica is going up or down
        going_up[top] = 0
        going_down[top] = 1
        going_up[bottom] = 1
        going_down[bottom] = 0

        up_count += going_up[this_trace]
        down_count += going_down[this_trace]
    bar.finish()

    f_up = up_count / (up_count + down_count + 0.1)
    f_down = down_count / (up_count + down_count + 0.1)
    return f_up, f_down


def get_progress_bar(label, n_steps):
    widgets = [
        "{}: ".format(label),
        progressbar.Percentage(),
        " ",
        progressbar.Bar(),
        " ",
        progressbar.ETA(),
    ]
    bar = progressbar.ProgressBar(maxval=n_steps, widgets=widgets)
    return bar


if __name__ == "__main__":
    main()
