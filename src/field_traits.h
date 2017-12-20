/*
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * Copyright (C) 2017 Nader Ganaba
 * Created by Nader on 19/12/2017.
 * field_traits.h
 */
#ifndef field_traits_h
#define field_traits_h

template< class Container >
struct field_traits{
    typedef Container container_type;
    typedef typename container_type::StateType StateType;
    typedef typename container_type::Iterator Iterator;
    
};
#endif /* field_traits_h */
