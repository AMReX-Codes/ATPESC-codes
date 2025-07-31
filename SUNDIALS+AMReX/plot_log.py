#!/usr/bin/env python3


def main():

    import argparse
    import matplotlib.pyplot as plt
    import matplotlib.ticker as tik

    from suntools import logs as sunlog

    parser = argparse.ArgumentParser(description="Plots")

    parser.add_argument("logfiles", type=str, nargs="+", help="Log files to plot")

    parser.add_argument("--logx", action="store_true", help="Use log scale for x-axis")

    parser.add_argument("--logy", action="store_true", help="Use log scale for y-axis")

    parser.add_argument("--labels", type=str, nargs="+", help="Plot labels")

    parser.add_argument(
        "--save",
        type=str,
        nargs="?",
        const="fig.pdf",
        default=None,
        metavar="FILE_NAME",
        help="Save figure to file",
    )

    # parse command line args
    args = parser.parse_args()

    _, ax = plt.subplots()
    colors = plt.get_cmap("tab10")

    for idx, lf in enumerate(args.logfiles):

        # parse log file
        log = sunlog.log_file_to_list(lf)

        # get all step sizes attempted
        _, times_s, vals_s = sunlog.get_history(log, "h")

        # get which steps failed because of too much error
        _, times_etf, vals_etf = sunlog.get_history(log, "h", "failed error test")

        # get which steps failed because the nonlinear solver failed
        _, times_sf, vals_sf = sunlog.get_history(log, "h", "failed solve")

        # plot log data
        if len(args.logfiles) == 1:
            s_color = "blue"
            etf_color = "red"
            sf_color = "darkorange"
        else:
            s_color = colors(idx)
            etf_color = s_color
            sf_color = s_color

        if args.labels:
            s_label = f"{args.labels[idx]} step attempt"
            etf_label = f"{args.labels[idx]} error test failed"
            sf_label = f"{args.labels[idx]} solver failed"
        else:
            s_label = "step attempt"
            etf_label = "error test failed"
            sf_label = "solver failed"

        # plot successful data
        ax.plot(times_s, vals_s, color=s_color, marker=".", label=s_label, zorder=0.1)

        # add failures as scatter plot
        ax.scatter(times_etf, vals_etf, color=etf_color, marker="x", label=etf_label, zorder=0.2)

        if len(vals_sf) > 0:
            ax.scatter(times_sf, vals_sf, color=sf_color, marker="d", label=sf_label, zorder=0.2)

    if args.logx:
        ax.set_xscale("log")
    if args.logy:
        ax.set_yscale("log")

    ax.set_xlabel("time")
    ax.set_ylabel("step size")
    ax.legend(loc="best")

    ax.grid(alpha=0.3, linestyle="--")

    if args.save:
        plt.savefig(args.save, bbox_inches="tight")
    else:
        plt.show()


# run the main routine
if __name__ == "__main__":
    import sys

    sys.exit(main())
