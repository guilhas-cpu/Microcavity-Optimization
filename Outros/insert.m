function NewCharArray = insert( CharArray, Position, WhatToInsert)
  NewCharArray = char( strcat( cellstr(CharArray(:,1:Position)), cellstr(WhatToInsert), cellstr(CharArray(:, Position+1:end)) ) );