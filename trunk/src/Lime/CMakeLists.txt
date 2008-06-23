SET(Sources_clime c-lime/lime_fseeko.cpp c-lime/lime_header.cpp c-lime/lime_reader.cpp c-lime/lime_utils.cpp c-lime/lime_writer.cpp)
SET(Sources_Data Data/Data_a.cpp Data/Data_b.cpp Data/Data.h)
SET(Sources_Reader Reader/Reader_a.cpp Reader/Reader_b.cpp Reader/Reader_c.cpp Reader/Reader.h)
ADD_LIBRARY(Lime ${Sources_clime} ${Sources_Data} ${Sources_Reader})