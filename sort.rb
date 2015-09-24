f = File.open(ARGV.shift, "r")
palindromes = []
while line = f.gets
  left, right, size, rate = line.split(";").map{|i| i.to_i}
  palindromes << [left, right, size, rate]
end

palindromes.sort! {|x, y| x[0] <=> y[0]}

palindromes.each do |p|
  puts "#{p[0]};#{p[1]};#{p[2]};#{p[3]}"
end
