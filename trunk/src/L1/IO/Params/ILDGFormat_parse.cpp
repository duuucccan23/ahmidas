#include "Params.ih"

void IO::ILDGFormat::parse(char *message)
{
  // We use the C tokenize capabilities to parse this string
  char *pch;
  pch = std::strtok(message, "<>");
  if (std::strncmp(pch, "?xml", 4))
    return;

  // We've removed the XML header, now we can set up a state machine to parse the file
  while ((pch = std::strtok(0, "<>")))
  {
    if (std::strncmp(pch, "ildgFormat", 10))
    {
      ildgFormat.assign(realFront(pch + 10), realLen(pch)); // Skip the tag itself, copy URL.
      continue;
    }

    while ((pch = std::strtok(0, "<>")))
    {
      if (!std::strncmp(pch, "/ildgFormat", 11))
        break; // We're done with the ildgFormat block altogether.
        
      if (!std::strncmp(pch, "version", 7))
      {
        while ((pch = std::strtok(0, "<>")))
        {
          if (!std::strncmp(pch, "/version", 8))
            break;
          version.assign(realFront(pch), realLen(pch));
        }
        continue;
      }
      if (!std::strncmp(pch, "field", 5))
      {
        while ((pch = std::strtok(0, "<>")))
        {
          if (!std::strncmp(pch, "/field", 6))
            break;
          field.assign(realFront(pch), realLen(pch));
        }
        continue;
      }
      if (!std::strncmp(pch, "precision", 9))
      {
        while ((pch = std::strtok(0, "<>")))
        {
          if (!std::strncmp(pch, "/precision", 10))
            break;
          precision = atoi(pch);
        }
        continue;
      }
  
      if (!std::strncmp(pch, "lx", 2))
      {
        while ((pch = std::strtok(0, "<>")))
        {
          if (!std::strncmp(pch, "/lx", 3))
            break;
          nx = atoi(pch);
        }
        continue;
      }
      if (!std::strncmp(pch, "ly", 2))
      {
        while ((pch = std::strtok(0, "<>")))
        {
          if (!std::strncmp(pch, "/ly", 3))
            break;
          ny = atoi(pch);
        }
        continue;
      }
      if (!std::strncmp(pch, "lz", 2))
      {
        while ((pch = std::strtok(0, "<>")))
        {
          if (!std::strncmp(pch, "/lz", 3))
            break;
          nz = atoi(pch);
        }
        continue;
      }
      if (!std::strncmp(pch, "lt", 2))
      {
        while ((pch = std::strtok(0, "<>")))
        {
          if (!std::strncmp(pch, "/lt", 3))
            break;
          nt = atoi(pch);
        }
        continue;
      }
    }
    break;
  }
}
