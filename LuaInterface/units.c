static int PhysicalConstantsFailed (lua_State *L)
{
  char *call;
  call = (char*)luaL_checkstring(L, -1);
  return luaL_error(L, "PhysicalConstants does not have the element: %s",call);
}

static int AllPhysicalConstantsToString (lua_State *L)
{
  const uint32_t nelements=sizeof(PhysicalConstants)/sizeof(PhysicalConstants[0]);
  uint32_t i,buffsize=20;
  char valuebuff[50];
  char *buff;
  for(i=0;i<nelements;i++)
  {
    buffsize+=strlen(PhysicalConstants[i].name)+strlen(PhysicalConstants[i].units)+strlen(PhysicalConstants[i].description)+50;
  }
  buff=malloc(buffsize*sizeof(char));
  strcpy(buff, "PhysicalConstants:\n\n");
  for(i=0;i<nelements;i++)
  {
    strcat(buff, "-- ");
    strcat(buff,PhysicalConstants[i].description);
    strcat(buff, "\n");
    sprintf(valuebuff, "%-8s",PhysicalConstants[i].name);
    strcat(buff, valuebuff);
    strcat(buff, "= ");
    sprintf(valuebuff, "%22.15E",PhysicalConstants[i].value);
    strcat(buff, valuebuff);
    strcat(buff, "  --[[ ");
    strcat(buff,PhysicalConstants[i].units);
    strcat(buff, " ]]--\n\n");
    
  }
  lua_pushstring(L, buff);
  free(buff);
  return 1;
}

static int PhysicalConstantsToString (lua_State *L)
{
  char *units, *description, *buff;
  double value;
  lua_pushstring(L, "value");
  lua_rawget(L, -2);
  value = luaL_checknumber(L, -1);
  lua_pop(L, 1);
  lua_pushstring(L, "units");
  lua_rawget(L, -2);
  units = (char*)luaL_checkstring(L, -1);
  lua_pop(L, 1);
  lua_pushstring(L, "description");
  lua_rawget(L, -2);
  description = (char*)luaL_checkstring(L, -1);
  lua_pop(L, 1);
  buff=malloc((strlen(description)+strlen(units)+50)*sizeof(char));
  sprintf(buff, "-- %s\n%22.15E --[[ %s ]]--",description,value,units);
  lua_pushstring(L, buff);
  free(buff);
  return 1;
}

static int RegisterPhysicalConstants(lua_State *L)
{
  const uint32_t nelements=sizeof(PhysicalConstants)/sizeof(PhysicalConstants[0]);
  uint32_t i;
  
  lua_createtable(L, 0, nelements);
  for(i=0;i<nelements;i++)
  {
    lua_pushstring(L, PhysicalConstants[i].name);
    lua_createtable(L, 0, 3);
    lua_pushstring(L, "value");
    lua_pushnumber(L, PhysicalConstants[i].value);
    lua_rawset(L, -3);
    lua_pushstring(L, "units");
    lua_pushstring(L, PhysicalConstants[i].units);
    lua_rawset(L, -3);
    lua_pushstring(L, "description");
    lua_pushstring(L, PhysicalConstants[i].description);
    lua_rawset(L, -3);
    
    lua_createtable(L, 0, 1);
    lua_pushcfunction(L, PhysicalConstantsToString);
    lua_setfield(L, -2, "__tostring");
    lua_setmetatable(L, -2);
    
    lua_settable(L, -3);
  }
  lua_createtable(L, 0, 2);
  lua_pushcfunction(L, PhysicalConstantsFailed);
  lua_setfield(L, -2, "__index");
  lua_pushcfunction(L, AllPhysicalConstantsToString);
  lua_setfield(L, -2, "__tostring");
  lua_setmetatable(L, -2);
  
  lua_setglobal(L, "PhysicalConstants");
  return 1;
}

static int EnergyUnitsNewIndex (lua_State *L)
{
  char *buff;
  buff=(char*)luaL_checkstring(L, -2);
  if(strcmp("Set",buff)==0)
  {
    buff=(char*)luaL_checkstring(L, -1);
    return RegisterEnergyUnits(L, buff);
  }
  else
  {
    return luaL_error(L, "You can change the standard Units in EnergyUnits by EnergyUnits.Set, no other options are available: %s",buff);
  }
}

static int EnergyUnitsIndex (lua_State *L)
{
  char *call;
  call = (char*)luaL_checkstring(L, -1);
  return luaL_error(L, "EnergyUnits does not have the element: %s",call);
}

static int AllEnergyUnitsToString (lua_State *L)
{
  uint32_t i,j;
  char buff[19];
  const uint32_t nelements=sizeof(EnergyUnits)/sizeof(EnergyUnits[0]);
  lua_pushstring(L, "Conversion table between energy units\n      ");
  for(i=0;i<nelements;i++)
  {
    sprintf(buff, "%-18s",EnergyUnits[i].namelong);
    lua_pushstring(L, buff);
  }
  lua_pushstring(L, "\n      ");
  lua_concat(L, nelements+2);
  for(i=0;i<nelements;i++)
  {
    sprintf(buff, "%-18s",EnergyUnits[i].nameshort);
    lua_pushstring(L, buff);
  }
  lua_pushstring(L, "\n");
  lua_concat(L, nelements+2);
  for(i=0;i<nelements;i++)
  {
    sprintf(buff, "%-6s",EnergyUnits[i].nameshort);
    lua_pushstring(L, buff);
    for(j=0;j<nelements;j++)
    {
      sprintf(buff, "%17.11E ",EnergyUnits[i].value/EnergyUnits[j].value);
      lua_pushstring(L, buff);
    }
    lua_pushstring(L, "\n");
    lua_concat(L, nelements+3);
  }
  return 1;
}

static int RegisterEnergyUnits(lua_State *L, char* EnUnits)
{
  const uint32_t nelements=sizeof(EnergyUnits)/sizeof(EnergyUnits[0]);
  uint32_t i,iEnUnits=0;
  double conversion=1;
  Bool found=false;
  i=0;
  while(!found)
  {
    if(i==nelements)
    {
      return luaL_error(L, "EnergyUnits tried to set standard unit: %s, but this unit does not excist",EnUnits);
    }
    if((strcmp(EnergyUnits[i].namelong, EnUnits)==0)||(strcmp(EnergyUnits[i].nameshort, EnUnits)==0))
    {
      found=True;
      conversion=1/EnergyUnits[i].value;
      iEnUnits=i;
    }
    i++;
  }
  
  lua_createtable(L, 0, 2*nelements);
  for(i=0;i<nelements;i++)
  {
    lua_pushstring(L, EnergyUnits[i].namelong);
    lua_createtable(L, 0, 3);
    lua_pushstring(L, "value");
    lua_pushnumber(L, conversion*EnergyUnits[i].value);
    lua_rawset(L, -3);
    lua_pushstring(L, "units");
    lua_pushstring(L, EnergyUnits[iEnUnits].nameshort);
    lua_rawset(L, -3);
    lua_pushstring(L, "description");
    lua_pushstring(L, EnergyUnits[i].description);
    lua_pushstring(L, " in units of ");
    lua_pushstring(L, EnergyUnits[iEnUnits].namelong);
    lua_concat(L, 3);
    lua_rawset(L, -3);
    lua_createtable(L, 0, 1);
    lua_pushcfunction(L, PhysicalConstantsToString);
    lua_setfield(L, -2, "__tostring");
    lua_setmetatable(L, -2);
    lua_settable(L, -3);
    
    lua_pushstring(L, EnergyUnits[i].nameshort);
    lua_createtable(L, 0, 3);
    lua_pushstring(L, "value");
    lua_pushnumber(L, conversion*EnergyUnits[i].value);
    lua_rawset(L, -3);
    lua_pushstring(L, "units");
    lua_pushstring(L, EnergyUnits[iEnUnits].nameshort);
    lua_rawset(L, -3);
    lua_pushstring(L, "description");
    lua_pushstring(L, EnergyUnits[i].description);
    lua_pushstring(L, " in units of ");
    lua_pushstring(L, EnergyUnits[iEnUnits].namelong);
    lua_concat(L, 3);
    lua_rawset(L, -3);
    lua_createtable(L, 0, 1);
    lua_pushcfunction(L, PhysicalConstantsToString);
    lua_setfield(L, -2, "__tostring");
    lua_setmetatable(L, -2);
    lua_settable(L, -3);
    
  }
  lua_createtable(L, 0, 3);
  lua_pushcfunction(L, EnergyUnitsNewIndex);
  lua_setfield(L, -2, "__newindex");
  lua_pushcfunction(L, EnergyUnitsIndex);
  lua_setfield(L, -2, "__index");
  lua_pushcfunction(L, AllEnergyUnitsToString);
  lua_setfield(L, -2, "__tostring");
  lua_setmetatable(L, -2);
  
  lua_setglobal(L, "EnergyUnits");
  return 1;
}

int luaopen_units(lua_State * L)
{
  RegisterPhysicalConstants(L);
  RegisterEnergyUnits(L, "eV");
  return 2;
}