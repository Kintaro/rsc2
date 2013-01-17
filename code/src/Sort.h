#ifndef __SORT_H__
#define __SORT_H__

class Sort
{
private:
	template<typename K, typename V>
	struct sort_helper 
	{
		static bool sort(const std::pair<K, V>& a, const std::pair<K, V>& b)  
		{
			return a.second < b.second;
		}
	};
public:
	template<typename K, typename V>
	static void partial_sort(std::vector<K>& a, std::vector<V>& b, const int from, const int to)
	{
		std::vector<std::pair<K, V>> temp;
		temp.resize(a.size());

		for (auto i = 0; i < a.size(); ++i)
			temp[i] = std::make_pair(a[i], b[i]);

		std::partial_sort(temp.begin() + 0, temp.begin() + to, temp.end(), sort_helper<K, V>::sort);

		for (auto i = 0; i < temp.size(); ++i)
		{
			a[i] = temp.first;
			b[i] = temp.second;
		}
	}
};

#endif