#!/usr/bin/env ruby
#  Copyright (c) 2013 University of Pennsylvania
#
#  Permission is hereby granted, free of charge, to any person obtaining a 
#  copy of this software and associated documentation files (the "Software"), 
#  to deal in the Software without restriction, including without limitation 
#  the rights to use, copy, modify, merge, publish, distribute, sublicense, 
#  and/or sell copies of the Software, and to permit persons to whom the 
#  Software is furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included in 
#  all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
#  DEALINGS IN THE SOFTWARE.

# for each locus compute the entropy of the read positions given an input of:
# locus   5p/3p   count

if ARGV.size < 1
  $stderr.puts "USAGE: #{$0} entropy_input"
  exit 1
end

infn = ARGV.shift

Locus = Struct.new(:name, :counts5p, :counts3p)

loci = Hash.new
File.open(infn).each_line do |line|
  name, which, count = line.chomp.split(/\t/)
  loci[name] ||= Locus.new(name, [], [])
  if which == "5p"
    loci[name].counts5p << count.to_i
  else
    loci[name].counts3p << count.to_i
  end
end

class Array
  def sum
    self.inject(0) { |a,x| a + x }
  end
end

puts "name\tpos_entropy5p\tpos_entropy3p"

loci.sort.each do |name, locus|
  count5p_sum = locus.counts5p.sum.to_f
  count3p_sum = locus.counts5p.sum.to_f
  probs5p = locus.counts5p.map { |x| x.to_f / count5p_sum }
  probs3p = locus.counts3p.map { |x| x.to_f / count3p_sum }
  entropy5p = -( probs5p.map { |prob| prob*Math.log(prob) }.sum )
  entropy3p = -( probs3p.map { |prob| prob*Math.log(prob) }.sum )
  puts "#{name}\t#{entropy5p}\t#{entropy3p}"
end
