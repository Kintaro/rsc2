#include "FileUtil.h"

/*-----------------------------------------------------------------------------------------------*/
bool FileUtil::read_option;
bool FileUtil::is_binary;
/*-----------------------------------------------------------------------------------------------*/
void FileUtil::open_read(const std::string& path, std::fstream& file, const boost::optional<bool>& binary)
{
	if (*binary)
		file.open(path.c_str(), std::ios::in | std::ios::binary);
	else
		file.open(path.c_str(), std::ios::in);
}
/*-----------------------------------------------------------------------------------------------*/
void FileUtil::open_write(const std::string& path, std::fstream& file, const boost::optional<bool>& binary)
{
	if (*binary)
		file.open(path.c_str(), std::ios::out | std::ios::binary);
	else
		file.open(path.c_str(), std::ios::out);
}
/*-----------------------------------------------------------------------------------------------*/
const void FileUtil::space(std::fstream& file)
{
	if (!Options::get_option_as<bool>("use-binary-files"))
		file << " ";
}
/*-----------------------------------------------------------------------------------------------*/
const void FileUtil::newline(std::fstream& file)
{
	if (!Options::get_option_as<bool>("use-binary-files"))
		file << std::endl;
}
/*-----------------------------------------------------------------------------------------------*/
