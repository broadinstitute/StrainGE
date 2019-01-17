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


import sys
import textwrap
import argparse
from abc import ABCMeta, abstractmethod


class Subcommand(metaclass=ABCMeta):
    """Represents a subcommand with its own argument parser arguments."""

    def register_arguments(self, subparser: argparse.ArgumentParser):
        """This function should register all required arguments for this
        subcommand."""

    @abstractmethod
    def __call__(self, *args, **kwargs):
        """When the subcommand is used on the command line, this function
        will be called."""


class SubcommandRegistry:
    def __init__(self, version=None, subcommands_title="", *args, **kwargs):
        self.parser = argparse.ArgumentParser(*args, **kwargs)
        self.parser.set_defaults(subcommand_func=None)
        self.subparsers = self.parser.add_subparsers(
            title=subcommands_title if subcommands_title else "Subcommands")

        if version:
            self.parser.add_argument(
                '--version', action='version', version=version)

    def register_subcommand(self, name: str, subcommand: Subcommand, **kwargs):
        # Use subcommand class level documentation also for documentation on
        # command line -h/--help
        if hasattr(subcommand.__class__, '__doc__'):
            subcommand_doc = subcommand.__class__.__doc__
            first_help_line = subcommand_doc.strip().split('\n\n')[0].strip()

            kwargs['help'] = first_help_line
            kwargs['description'] = textwrap.dedent(subcommand_doc)
            kwargs['formatter_class'] = argparse.RawDescriptionHelpFormatter

        subparser = self.subparsers.add_parser(name, **kwargs)

        # Initialize subcommand arguments
        subcommand.register_arguments(subparser)
        subparser.set_defaults(subcommand_func=subcommand)

    def run(self, parser_args: argparse.Namespace):
        args_dict = vars(parser_args)
        subcommand_func = args_dict.pop('subcommand_func')

        if subcommand_func:
            rc = subcommand_func(**args_dict)
        else:
            self.parser.print_help()
            rc = 1

        if rc is None:
            rc = 0

        sys.exit(rc)
