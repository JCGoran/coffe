#!/usr/bin/env python3

from __future__ import annotations

from argparse import ArgumentParser

import numpy as np

from coffe import Coffe


def parse_args():
    parser = ArgumentParser(
        prog="COFFE",
        description="The command line interface to the COFFE code",
    )

    parser.add_argument(
        "-s",
        "--settings",
        type=str,
        help="The settings file to use",
        required=True,
    )

    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="In which file to store the output",
        required=True,
    )

    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s 3.0.0",
    )

    parser.add_argument(
        "-k",
        "--kind",
        type=str,
        help="The kind of output to compute.",
        choices=["bg", "cf", "mp", "cov", "avg_cov", "avg_mp"],
        required=True,
    )

    args = parser.parse_args()

    return args


def compute_output(args):
    cosmo = Coffe.from_file(args.settings)

    savetxt_kwargs = {
        "fmt": "%.10e",
        "delimiter": "\t",
        "encoding": "utf-8",
    }

    if args.kind == "bg":
        cosmo._background_init()

        output = np.transpose(
            [
                [
                    z,
                    cosmo.scale_factor(z),
                    cosmo.hubble_rate(z),
                    cosmo.hubble_rate_conformal(z),
                    cosmo.hubble_rate_conformal_derivative(z),
                    cosmo.growth_factor(z),
                    cosmo.growth_rate(z),
                    cosmo.comoving_distance(z),
                ]
                for z in np.linspace(0, 15, cosmo.background_sampling)
            ]
        )

        np.savetxt(
            args.output,
            np.transpose(output),
            header="\t".join(
                [
                    "1)z",
                    "2)a",
                    "3)H[1/Mpc]",
                    "4)conformal_H[1/Mpc]",
                    "5)conformal_H_prime[1/Mpc^2]",
                    "6)D1",
                    "7)f",
                    "8)comoving_distance[Mpc]",
                ]
            ),
            **savetxt_kwargs,
        )

    elif args.kind == "cf":
        result = cosmo.compute_corrfunc_bulk()
        np.savetxt(
            args.output,
            np.array([[_.r, _.mu, _.z, _.value] for _ in result]),
            header="\t".join(["r[Mpc]", "mu", "z", "corrfunc_value"]),
            **savetxt_kwargs,
        )

    elif args.kind == "mp":
        result = cosmo.compute_multipoles_bulk()
        np.savetxt(
            args.output,
            np.array([[_.r, _.l, _.z, _.value] for _ in result]),
            header="\t".join(["r[Mpc]", "l", "z", "multipoles_value"]),
            fmt=["%.10e", "%d", "%.10e", "%.10e"],
            **{key: value for key, value in savetxt_kwargs.items() if key != "fmt"},
        )
    elif args.kind == "avg_mp":
        result = cosmo.compute_average_multipoles_bulk()
        np.savetxt(
            args.output,
            np.array([[_.r, _.l, _.z_min, _.z_max, _.value] for _ in result]),
            header="\t".join(["r[Mpc]", "l", "z_min", "z_max", "avg_multipoles_value"]),
            fmt=["%.10e", "%d", "%.10e", "%.10e", "%.10e"],
            **{key: value for key, value in savetxt_kwargs.items() if key != "fmt"},
        )

    elif args.kind == "cov":
        result = cosmo.compute_covariance_bulk()
        np.savetxt(
            args.output,
            np.array([[_.r1, _.r2, _.l1, _.l2, _.z, _.value] for _ in result]),
            header="\t".join(
                ["r1[Mpc]", "r2[Mpc]", "l1", "l2", "z", "covariance_multipoles_value"]
            ),
            fmt=["%.10e", "%.10e", "%d", "%d", "%.10e", "%.10e"],
            **{key: value for key, value in savetxt_kwargs.items() if key != "fmt"},
        )

    elif args.kind == "avg_cov":
        result = cosmo.compute_average_covariance_bulk()
        np.savetxt(
            args.output,
            np.array(
                [[_.r1, _.r2, _.l1, _.l2, _.z_min, _.z_max, _.value] for _ in result]
            ),
            header="\t".join(
                [
                    "r1[Mpc]",
                    "r2[Mpc]",
                    "l1",
                    "l2",
                    "z_min",
                    "z_max",
                    "covariance_avg_multipoles_value",
                ]
            ),
            fmt=["%.10e", "%.10e", "%d", "%d", "%.10e", "%.10e", "%.10e"],
            **{key: value for key, value in savetxt_kwargs.items() if key != "fmt"},
        )


def main():
    args = parse_args()

    compute_output(args)


if __name__ == "__main__":
    main()
