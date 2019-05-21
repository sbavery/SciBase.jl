export
    readSTKfile,
    STKParseTime


function readSTKfile(filename)
    f = open(filename)
    lines = readlines(f)
    close(f)

    colInds = STKInds(lines)

    dict = Dict()
    newDict = false
    currDict, currTitles = "", ""

    for i = 3:length(lines)
        if length(split(lines[i])) == 0 || strip(lines[i])[1] == '-'
            continue
        end
        if in(split(lines[i])[1], ["Min", "Max", "Total", "Mean", "Global", "No"])
            continue
        end

        if length(split(lines[i])) == 1
            dict[strip(lines[i])] = Dict()
            currDict = dict[strip(lines[i])]
            newDict = true
        else
            if newDict
                currTitles = STKParseTitles!(currDict, colInds, lines[i])
                newDict = false
            else
                if currDict == ""
                    dict[strip(lines[2])] = Dict()
                    currDict = dict[strip(lines[2])]
                    currTitles = STKParseTitles!(currDict, colInds, lines[i])
                else
                    if strip(lines[i][colInds[1]:colInds[2]-1]) == currTitles[1]
                        continue
                    else
                        if strip(lines[i][colInds[1]:colInds[2]-1])[end] == currTitles[1][end]
                            currTitles = STKParseTitles!(currDict, colInds, lines[i])
                        else
                            STKParseLine!(currDict, currTitles, colInds, lines[i])
                        end
                    end
                end
            end
        end
    end
    return dict
end


function STKInds(lines)
    line = ""
    for i = 3:length(lines)
        if length(split(lines[i])) > 1
            if split(lines[i])[1][1] == '-' && split(lines[i])[2][1] == '-'
                line = lines[i]
                break
            end
        end
    end

    colInds = Int64[]
    wasDash = false
    for j = 1:length(line)
        if line[j] == '-'
            if !wasDash
                wasDash = true
                push!(colInds, j)
            end
        else
            wasDash = false
        end
    end
    return [colInds; length(line)+1]
end


function STKParseTitles!(dict, colInds, line)
    keys = String[]
    for i = 1:length(colInds)-1
        key = strip(line[colInds[i]:colInds[i+1]-1])
        dict[key] = []
        push!(keys, key)
    end
    return keys
end


function STKParseLine!(dict, titles, colInds, line)
    for i = 1:length(titles)
        val = strip(line[colInds[i]:colInds[i+1]-1])
        if contains(titles[i], "(UTCG)")
            push!(dict[titles[i]], STKParseTime(val))
        elseif contains(titles[i], "sec") || contains(titles[i], "Percent")
            push!(dict[titles[i]], parse(Float64, val))
        elseif titles[i] == "Access"
            push!(dict[titles[i]], parse(Int64, val))
        else
            push!(dict[titles[i]], val)
        end
    end
end


function STKParseTime(s) # Takes in an STK time string e.g. "23 Dec 2015 23:50:33.123" and parses it into a Julia DateTime.
    s = strip(s)
    if length(s) == 23
        s = "0" * s
    elseif s == "N/A"
        return NaN
    end
    month = Dates.abbrenglish[lowercase(s[4:6])]
    time = DateTime(parse(Int, s[8:11]), month, parse(Int, s[1:2]), parse(Int, s[13:14]), parse(Int, s[16:17]), parse(Int, s[19:20]), parse(Int, s[22:24]))
end