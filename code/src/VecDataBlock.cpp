#include "VecDataBlock.h"

/*-----------------------------------------------------------------------------------------------*/
const size_t VecDataBlock::get_global_offset() const
{
	return this->global_offset;
}
/*-----------------------------------------------------------------------------------------------*/
const void VecDataBlock::set_global_offset(const size_t new_offset)
{
	this->global_offset = new_offset;
}
/*-----------------------------------------------------------------------------------------------*/
const size_t VecDataBlock::get_number_of_items() const
{
	return this->number_of_items;
}
/*-----------------------------------------------------------------------------------------------*/
const std::string VecDataBlock::get_filename_prefix() const
{
	return "foo";
}
/*-----------------------------------------------------------------------------------------------*/
