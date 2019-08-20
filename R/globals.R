# Copyright (c) 2019 The President and Fellows of Harvard College

if(getRversion() >= "2.15.1"){
  # To remove 'no visible binding for global variable ...' in data.table commands
  utils::globalVariables(c('t0', 'id', 'J', 'k'))
}
