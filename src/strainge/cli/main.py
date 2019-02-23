#  Copyright (c) 2016-2019, Broad Institute, Inc. All rights reserved.
#
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
#  * Neither the name Broad Institute, Inc. nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
#  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
#  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
#  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
#  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#  POSSIBILITY OF SUCH DAMAGE.
#

import logging
import argparse

import strainge
from strainge.cli.registry import SubcommandRegistry
from strainge.cli.kmersets import (KmerizeSubcommand, KmersimSubCommand,
                                   ClusterSubcommand, CreateDBSubcommand,
                                   PlotSubcommand, StatsSubcommand)
from strainge.cli.straingst import StrainGSTSubCommand
from strainge.cli.straingr import (CallSubcommand, ViewSubcommand,
                                   CompareSubCommand)

logger = logging.getLogger()


class StrainGECLI(SubcommandRegistry):
    """
    StrainGE main entry point class.

    Collects all available subcommands and builds a argument parser.
    """

    def __init__(self):
        desc = strainge.__doc__
        desc += "\n\nVersion: {}".format(strainge.__version__)

        super().__init__(description=desc, version=strainge.__version__,
                         formatter_class=argparse.RawDescriptionHelpFormatter)

        # Add arguments to control log level/verboseness
        self.parser.add_argument(
            '-v', '--verbose', action='count', default=0, required=False,
            help="Increase verbosity level, number of levels: 0, 1, 2"
        )

    def __call__(self, *args, **kwargs):
        """This is basically our main() function. Setup logging, determine
        which subcommand is called, and run the corresponding `Subcommand`
        instance."""

        args = self.parser.parse_args()

        # Setup logging
        logger.setLevel(logging.INFO)
        spaghetti_logger = logging.getLogger('spaghetti')
        spaghetti_logger.setLevel(logging.WARNING)

        formatter = logging.Formatter(
            "%(asctime)s - %(levelname)s:%(name)s:%(message)s")
        handler = logging.StreamHandler()
        handler.setLevel(logging.DEBUG)
        handler.setFormatter(formatter)
        logger.addHandler(handler)

        if args.verbose > 0:
            spaghetti_logger.setLevel(logging.INFO)

        if args.verbose > 1:
            logger.setLevel(logging.DEBUG)
            spaghetti_logger.setLevel(logging.DEBUG)

        self.run(args)


strainge_cli = StrainGECLI()

strainge_cli.register_subcommand('kmerize', subcommand=KmerizeSubcommand())
strainge_cli.register_subcommand('kmersim', subcommand=KmersimSubCommand())
strainge_cli.register_subcommand('cluster', subcommand=ClusterSubcommand())
strainge_cli.register_subcommand('createdb', subcommand=CreateDBSubcommand())

strainge_cli.register_subcommand('search', subcommand=StrainGSTSubCommand())

strainge_cli.register_subcommand('call', subcommand=CallSubcommand())
strainge_cli.register_subcommand('view', subcommand=ViewSubcommand())
strainge_cli.register_subcommand('compare', subcommand=CompareSubCommand())

strainge_cli.register_subcommand('stats', subcommand=StatsSubcommand())
strainge_cli.register_subcommand('plot', subcommand=PlotSubcommand())
