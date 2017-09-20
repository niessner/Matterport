function loadMatterport(file_name, root_path)
  print(string.format("Loading data from %s", file_name))
  local output = {}
  local file = io.open(file_name, "r")
  for line in file:lines() do
    line = root_path .. line 

    local data = {}
    data.color = line
    baseline = line
    baseline = baseline:gsub("_i", "_d")
    baseline = baseline:gsub("undistorted_color_dmages", "undistorted_normal_images")
    
    data.nx = string.gsub(baseline, '.jpg', '_nx.png')
    data.ny = string.gsub(baseline, '.jpg', '_ny.png')
    data.nz = string.gsub(baseline, '.jpg', '_nz.png')

    
    data.name = string.gsub(baseline, '.jpg', '')

    table.insert(output, data)
  end
  print(string.format("%d data loaded.",#output))
  return output
end
