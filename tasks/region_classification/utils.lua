
-- Check if a file exists
function fileExists(file)
    local f = io.open(file, "rb")
    if f then f:close() end
    return f ~= nil
end

-- Get all lines from a file, returns an empty list/table if the file does not exist
function getLinesFromFile(file)
    if not fileExists(file) then return {} end
    lines = {}
    for line in io.lines(file) do 
        lines[#lines + 1] = line
    end
    return lines
end

-- Shuffle table
function shuffleTable(t,n)
    while n > 2 do -- only run if the table has more than 1 element
        local k = math.random(n) -- get a random number
        t[n], t[k] = t[k], t[n]
        n = n - 1
    end
    return t
end

-- Lookup filenames in directory (with search query string)
function scanDir(directory,query)
    local i, t, popen = 0, {}, io.popen
    local pfile = popen('ls -a "'..directory..'"')
    for filename in pfile:lines() do
        if filename == '.' or  filename == '..' then 
        else
            if query then
                if string.find(filename,query) then
                    i = i+1
                    t[i] = filename
                end
            else
                print(filename)
                i = i+1
                t[i] = filename
            end
        end
    end
    pfile:close()
    return t
end