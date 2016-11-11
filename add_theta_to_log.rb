#!/usr/bin/env ruby

# Michael Matschiner, 2016-11-10
#
# This script reads .trees and .log files produced with SNAPP, extracts the theta
# estimates from the .trees file and adds these to the .log file. It is assumed that
# the XML input file for SNAPP was produced with the script snapp_prep.rb. In such XML
# files, all theta values are linked. When the species tree has many branches, this
# would lead to a lot of redundant output (and potentially large file size) in the log
# files if all theta estimates are written to these. Therefore, theta estimates are only 
# written to the .trees files. To re-insert a single column with theta estimates to the
# log file subsequent to the analysis, this script can be used.
# Instead of overwriting the input log file, a new log file is written.
#
# To learn about the available options of this script, run with '-h' (or without any options):
# ruby snapp_prep.rb -h
#
# For questions and bug fixes, email to michaelmatschiner@mac.com

# Load required libraries.
require 'optparse'

# Define default options.
options = {}
options[:log] = "snapp.log"
options[:trees] = "snapp.trees"
options[:generation_time] = 1.0
options[:out] = options[:log].chomp(".log") + "_w_theta.log"

# Get the command line options.
ARGV << '-h' if ARGV.empty?
opt_parser = OptionParser.new do |opt|
	opt.banner = "Usage: ruby #{$0} [OPTIONS]"
	opt.separator  ""
	opt.separator  "Example"
	opt.separator  "ruby #{$0} -l #{options[:log]} -t #{options[:trees]} -g #{options[:generation_time]} -o #{options[:out]}"
	opt.separator  ""
	opt.separator  "Options"
	opt.on("-l","--log FILENAME","Name of .log output file of SNAPP analysis (default: #{options[:log]}).") {|l| options[:log] = l}
	opt.on("-t","--trees FILENAME","Name of .trees output file of SNAPP analysis (default: #{options[:trees]}).") {|t| options[:trees] = t}
	opt.on("-g","--gentime NUMBER",Float,"Generation time (required to calculate population size) (default: #{options[:generation_time]}).") {|g| options[:generation_time] = g}
	opt.on("-o","--out FILENAME","Name of  new .log output file (default: #{options[:out]}).") {|o| options[:out] = o}
	opt.on("-h","--help","Print this help text.") {
		puts opt_parser
		exit(0)
	}
	opt.separator  ""
end
opt_parser.parse!

# Read the .log input file.
log_file =  File.open(options[:log])
log_lines = log_file.readlines
log_header_line = log_lines[0]
log_state_lines = log_lines[1..-1]
log_state_numbers = []
log_clock_rates = []
log_state_lines.each do |l|
	log_state_numbers << l.split[0]
	log_clock_rates << l.split[-1].to_f
end

# Read the .trees input file.
trees_file = File.open(options[:trees])
trees_lines = trees_file.readlines
tree_states = []
thetas = []
trees_lines.each do |l|
	if l.include?("STATE_")
		l.match(/STATE_(\d+)\s*=/)
		tree_state = $1
		if l.include?("theta")
			l.match(/theta=([0-9\.\-eE]+)/)
		elsif l.include?("null")
			l.match(/null=([0-9\.\-eE]+)/)
		end
		theta = $1.to_f
		if tree_state != nil and theta != nil
			tree_states << tree_state
			thetas << theta
		end
	end
end

# Add theta estimates to .log lines.
out_string = "#{log_header_line.strip}\ttheta\tpopulation_size\t\n"
log_state_numbers.size.times do |x|
	if tree_states.include?(log_state_numbers[x])
		theta = thetas[tree_states.index(log_state_numbers[x])]
		mutation_rate = log_clock_rates[x]/(1000000.0/options[:generation_time])		
		pop_size = theta/(4.0*mutation_rate)
		out_string << "#{log_state_lines[x].strip}\t#{theta}\t#{pop_size}\t\n"
	else
		out_string << "#{log_state_lines[x].strip}\tNA\tNA\t\n"
	end
end

# Write the SNAPP input file.
out_file = File.open(options[:out],"w")
out_file.write(out_string)
