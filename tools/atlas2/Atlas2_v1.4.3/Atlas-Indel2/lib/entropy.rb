#!/usr/bin/ruby

def entropy(a_str, window_size)
	len=a_str.length
	pattern = Array.new(2*len, 0)
	return 0.0 if(len < window_size)
	n = 0
	i=0
	while i+window_size < len
		found = false
		p=a_str[i...a_str.length]
		j=0
		while j < n
			match = true
			q = a_str[(pattern[j*2])...a_str.length]
			for k in 0...window_size
				if(p[k] != q[k])
					match=false
					break
				end
			end
			if(match)
				found=true
				pattern[j*2+1]+=1
				break
			end
			j+=1
		end
		if(!found)
			pattern[n*2]=i
			pattern[n*2+1]=1
			n += 1
		end
		i += 1
	end
	e=0
	s=1.0/(len-window_size+1)
	for i  in 0...n
		f=s*pattern[i*2+1]
		e-=f*Math.log(f)
	end
	return e
end
