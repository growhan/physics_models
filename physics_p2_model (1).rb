require "csv"
class PhysicsModel
  M1 = 100000
  M_PLANET = ( 1.1944 * (10 ** 25))
  G_CONST = 6.67408 * ( 10 ** (-11) )
  F_PROPULSION = 20000
  TIME_STEP = 1
  def initialize
    @variables
    @data_hash = { 0 => { r: ( 4 * (10 ** 6)), r_theta: 90, f_g: nil, f_g_theta: nil, f_p: F_PROPULSION, f_p_theta: nil, f_net_x: nil,
                        f_net_y: nil, f_net: nil, f_net_theta: nil, a: nil, v: 0, delta_x_prior: 0   }}
  end
  attr_accessor :data_hash
  def build_data
    current_time = 0
    while current_time < 100000
      puts "TIME => #{current_time}"
      build_row(current_time)
      print @data_hash[current_time][:r]
      print ", "
      print @data_hash[current_time][:r_theta]
      print " || "
      print @data_hash[current_time][:f_net]
      print ", "
      print @data_hash[current_time][:f_net_theta]
      print @data_hash[current_time][:r_theta]
      puts
      build_next_row
      current_time += TIME_STEP
      break if @data_hash[current_time - TIME_STEP][:r_theta] <= 0 
    end
  end
  def write_to_csv
    file_name = "physics_data_m1_" + M1.to_s + ".csv"
    CSV.open(file_name, "wb") do |csv|
      csv << [ "time ( second )", "x_coord", "y_coord", "com_x_coord", "com_y_coord", "f_net", "f_net direction", "acceleration", "velocity", "KE", "PE"]
      @data_hash.each do |key, value|
        x_coord = value[:r] * cos(value[:r_theta])
        y_coord = value[:r] * sin(value[:r_theta])
        com_x_coord = (x_coord * M1) / ( M1 + M_PLANET)
        com_y_coord = (y_coord * M1) / ( M1 + M_PLANET)
        ke = 0.5 * M1 * (( value[:v]) ** 2)
        pe = M1 * 9.8 * value[:r]
        print [key, x_coord, y_coord, com_x_coord, com_y_coord, value[:f_net], value[:f_net_theta], value[:a], value[:v], ke, pe ]
        print "\n"
        csv << [key, x_coord, y_coord, com_x_coord, com_y_coord, value[:f_net], value[:f_net_theta], value[:a], value[:v], ke, pe ]
      end
    end
  end
  def build_row(time)
    row = @data_hash[time]
    row[:f_g] = get_f_g(row[:r])
    row[:f_g_theta] = row[:r_theta] - 180
    row[:f_p_theta] = 90 - row[:r_theta]
    row[:f_net_x] = ( row[:f_p] * cos(row[:f_p_theta])) - ( row[:f_g] * cos(row[:f_g_theta]))
    row[:f_net_y] = ( row[:f_p] * sin(row[:f_p_theta])) - ( row[:f_g] * sin(row[:f_g_theta]))
    row[:f_net] = pythag_theor(row[:f_net_x], row[:f_net_y])
    row[:f_net_theta] = -1 * tan_inv(row[:f_net_x], row[:f_net_y])
    row[:a] = row[:f_net] / M1
  end
  def build_next_row
    key_value = @data_hash.keys().last + TIME_STEP
    row_before = @data_hash[key_value - TIME_STEP] 
    @data_hash[key_value] = {}
    row = @data_hash[key_value]
    row[:v] = row_before[:v] + (row_before[:a] * TIME_STEP)
    # delta_x_prior is magnitude of displacement vector
    row[:delta_x_prior] = (row[:v] * TIME_STEP) / 2
    row_x = (row_before[:r] * cos(row_before[:r_theta])) + (row[:delta_x_prior]* cos(row_before[:f_net_theta]))
    row_y = (row_before[:r] * sin(row_before[:r_theta])) + (row[:delta_x_prior]* sin(row_before[:f_net_theta]))
    row[:r] = pythag_theor(row_x, row_y)
    row[:r_theta] = tan_inv(row_x, row_y)
    row[:f_p] = F_PROPULSION
  end
  def get_f_g(r)
    ( G_CONST * M1 * M_PLANET ) / ( r ** 2)
  end
  def sin(x)
    Math.sin(x * (Math::PI / 180)).round(6)
  end
  def cos(x)
    Math.cos(x * (Math::PI / 180)).round(6)
  end
  def tan_inv(len1, len2)
    ((Math.atan(len2/len1)) * ( 180 / Math::PI )).round(6)
  end
  def pythag_theor(len1, len2)
    ((( len1 ** 2 ) + ( len2 ** 2)) ** 0.5).round(6)
  end
end

phy = PhysicsModel.new
phy.build_data
phy.write_to_csv