#!/usr/bin/env python3

from pycallgraph import PyCallGraph
from pycallgraph.output import GraphvizOutput
from pycallgraph import Config
from pycallgraph import GlobbingFilter
import key_validation
import config_parser

##################
# Parser section
##################

parser_graphviz = GraphvizOutput(output_file='../../media/CallGraphDiagrams/parser_call_graph.png')

parser_config = Config()
parser_config.trace_filter = GlobbingFilter(include=['config_parser.*', "util*"], exclude=['*listcomp*', "*genexpr*"])

with PyCallGraph(output=parser_graphviz, config=parser_config):
    config_parser.main(["-i", "config_parser_test_resources/test_input.config",
                       "--jsonTemplate", "config_parser_test_resources/test_output_template.json",
                       "-o", "config_parser_test_resources/test_output.json"
                        ]
                       )

#####################
# Validator section
#####################

validator_graphviz = GraphvizOutput(output_file='../../media/CallGraphDiagrams/validator_call_graph.png')

validator_config = Config()
validator_config.trace_filter = GlobbingFilter(include=['key_validation.*', "util*"])


with PyCallGraph(output=validator_graphviz, config=validator_config):
    key_validation.main(["-i", "validator_test_resources/test_config.json",
                         "--KeyTypeFile", "validator_test_resources/test_key_types.json"
                         ]
                        )
