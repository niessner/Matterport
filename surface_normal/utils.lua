function loadData(file_name, config)
	print(string.format("Loading data from %s", file_name))
	local output = {}
	local file = io.open(file_name, "r")
	for line in file:lines() do
		-- print(line)
	    local data = {}

	    -- define data here, will be used in BatchIterator
	    -- data.input_file  = line .. '_mlt.png'
      data.input_file  = line .. '_color.png'
	    data.input_valid = line .. '_valid.png'
	    data.output_file = line .. '_norm_camera.png'
	    data.name       = line
	    -- end
	    
	    table.insert(output, data)
	end
	print(string.format("%d data loaded.",#output))
  return output
end

function round(num, idp)
  local mult = 10^(idp or 0)
  return math.floor(num * mult + 0.5) / mult
end

function split(str, pat)
   local t = {}  -- NOTE: use {n = 0} in Lua-5.0
   local fpat = "(.-)" .. pat
   local last_end = 1
   local s, e, cap = str:find(fpat, 1)
   while s do
      if s ~= 1 or cap ~= "" then
	 table.insert(t,cap)
      end
      last_end = e+1
      s, e, cap = str:find(fpat, last_end)
   end
   if last_end <= #str then
      cap = str:sub(last_end)
      table.insert(t, cap)
   end
   return t
end

function file_exists(name)
   local f=io.open(name,"r")
   if f~=nil then io.close(f) return true else return false end
end

function bool2num(var)
  return var and 1 or 0
end

function read_weight(file_name)
    print(string.format("Loading weight from %s", file_name))
    local file = io.open(file_name, "r")
    a = file:read("*all")
    b = string.split(a, '\n')
    w = torch.Tensor(table.getn(b))
    for i = 1,table.getn(b),1 do
      w[i] = b[i]
    end
    return w
end