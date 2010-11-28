#include <L1/Tool/IO.h>

inline char *realFront(char *string)
{
  return string + std::strspn(string, " \t");
}

inline size_t realLen(char *string)
{
  return std::strcspn(realFront(string), " \t");
}

/***/

Tool::IO::ScidacInfo Tool::IO::parseScidacInfo(char *scidacCStr)
{
  ScidacInfo result;

  // We use the C tokenize capabilities to parse this string
  char *pch;
  pch = std::strtok(scidacCStr, "<>");
  if (std::strncmp(pch, "?xml", 4))
    return ScidacInfo();

  // We've removed the XML header, now we can set up a state machine to parse the file
  while ((pch = std::strtok(0, "<>")))
  {
    if (std::strncmp(pch, "etmcFormat", 10))
      continue; // Unknown info string, just go on

    while ((pch = std::strtok(0, "<>")))
    {
      if (!std::strncmp(pch, "/etmcFormat", 11))
        break;
      if (!std::strncmp(pch, "field", 5))
      {
        while ((pch = std::strtok(0, "<>")))
        {
          if (!std::strncmp(pch, "/field", 6))
            break;
          result.field.assign(realFront(pch), realLen(pch));
        }
        continue;
      }
      if (!std::strncmp(pch, "precision", 9))
      {
        while ((pch = std::strtok(0, "<>")))
        {
          if (!std::strncmp(pch, "/precision", 10))
            break;
          result.precision.assign(realFront(pch), realLen(pch));
        }
        continue;
      }
      if (!std::strncmp(pch, "flavours", 8))
      {
        while ((pch = std::strtok(0, "<>")))
        {
          if (!std::strncmp(pch, "/flavours", 9))
            break;
          result.flavours.assign(realFront(pch), realLen(pch));
        }
        continue;
      }
      if (!std::strncmp(pch, "lx", 2))
      {
        while ((pch = std::strtok(0, "<>")))
        {
          if (!std::strncmp(pch, "/lx", 3))
            break;
          result.dims[Base::idx_X] = atoi(realFront(pch));
        }
        continue;
      }
      if (!std::strncmp(pch, "ly", 2))
      {
        while ((pch = std::strtok(0, "<>")))
        {
          if (!std::strncmp(pch, "/ly", 3))
            break;
          result.dims[Base::idx_Y] = atoi(realFront(pch));
        }
        continue;
      }
      if (!std::strncmp(pch, "lz", 2))
      {
        while ((pch = std::strtok(0, "<>")))
        {
          if (!std::strncmp(pch, "/lz", 3))
            break;
          result.dims[Base::idx_Z] = atoi(realFront(pch));
        }
        continue;
      }
      if (!std::strncmp(pch, "lt", 2))
      {
        while ((pch = std::strtok(0, "<>")))
        {
          if (!std::strncmp(pch, "/lt", 3))
            break;
          result.dims[Base::idx_T] = atoi(realFront(pch));
        }
        continue;
      }
    }
    break;
  }
  return result;
}
