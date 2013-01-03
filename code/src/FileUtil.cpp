#include "FileUtil.h"

/*-----------------------------------------------------------------------------------------------*/
bool FileUtil::read_option;
bool FileUtil::is_binary;
/*-----------------------------------------------------------------------------------------------*/
const void FileUtil::space(std::fstream& file)
{
	if (!read_option)
	{
		is_binary = Options::get_option_as<bool>("use-binary-files");
		read_option = true;
	}

	if (!is_binary)
		file << " ";
}
/*-----------------------------------------------------------------------------------------------*/
const void FileUtil::newline(std::fstream& file)
{
	if (!read_option)
	{
		is_binary = Options::get_option_as<bool>("use-binary-files");
		read_option = true;
	}

	if (!is_binary)
		file << std::endl;
}
/*-----------------------------------------------------------------------------------------------*/
