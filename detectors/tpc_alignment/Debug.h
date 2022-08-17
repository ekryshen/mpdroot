#ifndef DEBUG_HH
#define DEBUG_HH

namespace TpcAlignment {

	class Debug
	{
	private:
		bool vDebugMode{ false };

	public:
		Debug(bool mode);
		~Debug();
		/**
		* @brief Print debug message if allowed.
		*
		* @param _Format
		* @param ... - arguments for printf function
		*/
		void Print(char const *const _Format, ...) const;
	};
}

#endif