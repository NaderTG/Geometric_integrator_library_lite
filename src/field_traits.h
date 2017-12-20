//
//  field_traits.h
//  vorticity_lite
//
//  Created by Nader on 20/12/2017.
//  Copyright © 2017 Nader Ganaba. All rights reserved.
//

#ifndef field_traits_h
#define field_traits_h

template< class Container >
struct field_traits{
    typedef Container container_type;
    typedef typename container_type::StateType StateType;
    typedef typename container_type::Iterator Iterator;
    
};
#endif /* field_traits_h */
