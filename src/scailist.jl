module scailist
export Interval, overlap, ScAIList, find

struct Interval{T}
	start::Int
	stop::Int
	val::T
end

function overlap(self::Interval, other::Interval)::Bool
	self.start < other.stop && self.stop > other.start
end

struct ScAIList{T}
	intervals::Vector{Interval{T}}
	num_comps::Int
	comp_lens::Vector{Int}
	comp_idxs::Vector{Int}
	max_ends::Vector{Int}
end

function ScAIList{T}(input_intervals::Vector{Interval{T}}, min_cov_len::Int = 20)::ScAIList{T} where T
	max_comps = floor(Int, log2(float(length(input_intervals)))) + 1
	min_cov = div(min_cov_len, 2)

	# A compoenet is essentially a sublist
	num_comps = nothing # number of components
	comp_lens = Vector{Int}() # lengths of each component
	comp_idxs = Vector{Int}() # start positions of each component
	max_ends = Vector{Int}() # list of the max end positions

	# min_comp_len = max(64, min_cov_len) # min length of a component
	min_comp_len = 2

	input_len = length(input_intervals)
	sort!(input_intervals, by  = iv -> iv.start)
	
	decomposed = Vector{Interval{T}}()

	if input_len <= min_comp_len
		num_comps = 1
		push!(comp_lens, input_len)
		push!(comp_idxs, 1)
		append!(decomposed, input_intervals)
	else
		curr_comp = 0
		while curr_comp < max_comps && input_len - length(decomposed) > min_comp_len
			list1 = Vector{Interval{T}}()
			list2 = Vector{Interval{T}}()
			for i in eachindex(input_intervals)
				interval = input_intervals[i]
				j = 2
				cov = 1
				while j < min_comp_len && cov < min_cov && j + i < length(input_intervals)
					if input_intervals[j+i].stop >= interval.stop
						cov += 1
					end
					j += 1
				end
				if cov < min_cov
					push!(list1, interval)
				else
					push!(list2, interval)
				end
			end

			# Add the component info to ScAIList
			push!(comp_idxs, length(decomposed) + 1)
			push!(comp_lens, length(list1))
			curr_comp += 1
			
			if length(list2) <= min_comp_len || curr_comp == max_comps - 2
				if !isempty(list2)
					append!(decomposed, list1)
					push!(comp_idxs, length(decomposed) + 1)
					push!(comp_lens, length(list2))
					append!(decomposed, list2)
					curr_comp += 1
				end
			else
				append!(decompsed, list1)
				input_intervals = list2
			end

		end
		num_comps = curr_comp
	end
	
	# Augment with maxend
	for j in 1:num_comps
		comp_start = comp_idxs[j]
		comp_end = comp_start + comp_lens[j] - 1
		max_end = decomposed[comp_start].stop
		push!(max_ends, max_end)
		for iv in @view decomposed[comp_start + 1:comp_end]
			if iv.stop > max_end
				max_end = iv.stop
			end
			push!(max_ends, max_end)
		end
	end
	
	return ScAIList{T}(decomposed, num_comps,comp_lens, comp_idxs, max_ends)
end

# Binary search to find the right most index where intervals.start < query.stop
function upper_bound{T}(stop::Int, intervals::Vector{Interval{T}})::Union{Nothing, Int} where T
	right = length(intervals) + 1
	left = 1

	if intervals[right - 1].start < stop
		return right - 1
	elseif intervals[left].start >= stop
		return nothing
	end

	while right > 1
		half = div(right, 2)
		other_half = right - half
		probe = left + half
		other_left = left + other_half
		v = intervals[probe]
		right = half
		left = v.start < stop ? other_left : left
	end

	# Guarded at the top from ending on either extreme
	if intervals[left].start >= stop
		return left - 1
	else
		return left
	end
end

function find(self::ScAIList, start::Int, stop::Int)
	return Find(self, stop, start)
end

struct Find{T}
	inner::ScAIList{T}
	stop::Int
	start::Int
end

function Base.iterate(it::Find, state=(1, 1, true, false))
	offset, comp_num, find_offset, breaknow = state
	
	while comp_num <= it.inner.num_comps
		comp_start = it.inner.comp_idxs[comp_num]
		comp_end = comp_start + it.inner.comp_lens[comp_num]
		if it.inner.comp_lens[comp_num] > 15
			println("In sublist section")
			if find_offset
				tmp_offset = upper_bound(it.stop, @view it.inner.intervals[comp_start:comp_end])
				if isnothing(tmp_offset)
					comp_num += 1
					find_offset = true
					continue
				else
					offset = tmp_offset
				end
				offset += comp_start
				find_offset = false
			end

			while offset >= comp_start && it.inner.max_ends[offset] > start && !breaknow
				interval = it.intervals[offset]
				offset -= 1 # Do I need to check going under 0?
				if interval.stop > it.start
					return (interval, (offset, comp_num, find_offset, breaknow))
				end
			end
		else
			while offset < comp_end
				interval = it.inner.intervals[offset]
				offset += 1
				if interval.start < it.stop && interval.stop > it.start
					return (interval, (offset, comp_num, find_offset, breaknow))
				end
			end
		end
		breaknow = false
		find_offset = true
		comp_num += 1
	end
	return nothing
end

# For now just return them in non-sorted order
# function Base.iterate(it::ScAIList, state=1)
	# if state < length(it.intervals)
# end

end # module
