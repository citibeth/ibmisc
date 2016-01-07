
/** Acts like a std::function that is made up of many sub-std::fns */
template<class FunctionT>
class CompositeFn {
	std::vector<FunctionT> _fns;
public:
	void operator+=(FunctionT const &fn)
		{ _fns.push_back(fn); }

	void operator()
	{
		for (auto ii=_fns.begin(); ii != _fns.end(); ++ii) (*ii)();
		_fns.clear();
	}
};