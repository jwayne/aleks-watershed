/*
Aleks's master's thesis watershed algorithm with the T_h, T_l, T_e thresholds.
Further merging of oversegmented regions is performed depending on the `type`
input.

Algorithm:
    (1) CC with threshold T_h.
    (2) Simple watershed (Cousty 2009) considering only edges > T_l
    (3) Merge oversegmented regions according to `type`..

Type:
    0: Merge all pairs of segments with max edge weight between them > T_e (0.1)
    1: Merge all pairs of segments with at least one part having size <=
       limit_fn_avg( avg edge weight between them )
    2: Merge all pairs of segments with at least one part having size <=
       limit_fn_bup( max edge weight between them )
Note that for types 1-2, merges of >2 segments into the same resulting label
are made with sizes updated dynamically.  This is quick since computing size is O(1).

Code implementing the T_s threshold discussed in Aleks's master's thesis is
included, but is commented out, having been replaced by types 1-2.

Code by Aleksander Zlateski
*/

#include <boost/multi_array.hpp>
#include <boost/multi_array/types.hpp>
#include <memory>
#include <type_traits>

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstddef>
#include <cstdint>
#include <queue>
#include <vector>
#include <algorithm>
#include <tuple>
#include <map>
#include <list>
#include <vector>
#include <chrono>

#include <zi/disjoint_sets/disjoint_sets.hpp>


template < typename T > struct watershed_traits;

template <> struct watershed_traits<uint32_t>
{
    static const uint32_t high_bit = 0x80000000;
    static const uint32_t mask     = 0x7FFFFFFF;
    static const uint32_t visited  = 0x00001000;
    static const uint32_t dir_mask = 0x0000007F;
};

template <> struct watershed_traits<uint64_t>
{
    static const uint64_t high_bit = 0x8000000000000000LL;
    static const uint64_t mask     = 0x7FFFFFFFFFFFFFFFLL;
    static const uint64_t visited  = 0x0000000000001000LL;
    static const uint64_t dir_mask = 0x000000000000007FLL;
};

template < typename T >
using volume = boost::multi_array<T,3>;

template < typename T >
using affinity_graph = boost::multi_array<T,4>;

template < typename T >
using volume_ptr = std::shared_ptr<volume<T>>;

template < typename T >
using affinity_graph_ptr = std::shared_ptr<affinity_graph<T>>;


template < typename T >
inline bool read_from_file( const std::string& fname, T* data, std::size_t n )
{
    FILE* f = std::fopen(fname.c_str(), "rbXS");

    if ( !f )
    {
        return false;
    }

    std::size_t nread = std::fread(data, sizeof(T), n, f);

    std::fclose(f);

    return nread == n;
}

template < typename T >
inline affinity_graph_ptr<T>
read_affinity_graph_from_file( const std::string& fname,
                               std::size_t xsize,
                               std::size_t ysize,
                               std::size_t zsize )
{
    affinity_graph_ptr<T> aff(new affinity_graph<T>
                              (boost::extents[xsize][ysize][zsize][3],
                               boost::fortran_storage_order()));

    if ( !read_from_file(fname, aff->data(), xsize*ysize*zsize*3) )
    {
        throw 0;
    }

    return aff;
}


template < typename T >
inline bool
write_volume_to_file( const std::string& fname,
                      const volume_ptr<T>& vol )
{
    /* Write labels to a file. */

    std::ofstream f(fname.c_str(), (std::ios::out | std::ios::binary) );

    if ( !f )
    {
        return false;
    }

    f.write( reinterpret_cast<char*>(vol->data()),
             vol->shape()[0] * vol->shape()[1] * vol->shape()[2] * sizeof(T));

    return true;
}

template< typename ID, typename F >
using region_graph = std::vector<std::tuple<F,ID,ID>>;

template< typename ID, typename F >
using region_graph_ptr = std::shared_ptr<region_graph<ID,F>>;


template< typename ID, typename F, typename L, typename H >
inline std::pair<volume_ptr<ID>, ID>
simple_watershed( const affinity_graph_ptr<F>& aff_ptr,
                  const L& lowv,
                  const H& highv,
                  std::vector<std::size_t>& counts )
{
    /*
    Perform watershed as defined by Algorithm 1 in Cousty 2009.

    Input:
      aff_ptr = affinity graph of edge weights
      lowv = T_l
      highv = T_h
      counts = empty.  will fill with size of each segment
    Output:
      pair (vertex labels--z*y*x float, # labels)
      Note that `counts` is also updated.
    */

    typedef F float_type;
    typedef watershed_traits<ID> traits;

    float_type low  = static_cast<float_type>(lowv);
    float_type high = static_cast<float_type>(highv);

    std::ptrdiff_t xdim = aff_ptr->shape()[0];
    std::ptrdiff_t ydim = aff_ptr->shape()[1];
    std::ptrdiff_t zdim = aff_ptr->shape()[2];

    std::ptrdiff_t size = xdim * ydim * zdim;

    volume_ptr<ID> seg_ptr( new volume<ID>(boost::extents[xdim][ydim][zdim],
                                           boost::fortran_storage_order()));
    affinity_graph<F>& aff = *aff_ptr;
    volume<ID>&        seg = *seg_ptr;

    ID* seg_raw = seg_ptr->data();

    // for each vertex, find which direction the stream is going
    for ( std::ptrdiff_t z = 0; z < zdim; ++z )
        for ( std::ptrdiff_t y = 0; y < ydim; ++y )
            for ( std::ptrdiff_t x = 0; x < xdim; ++x )
            {
                ID& id = seg[x][y][z] = 0;

                F negx = (x>0) ? aff[x][y][z][0] : low;
                F negy = (y>0) ? aff[x][y][z][1] : low;
                F negz = (z>0) ? aff[x][y][z][2] : low;
                F posx = (x<(xdim-1)) ? aff[x+1][y][z][0] : low;
                F posy = (y<(ydim-1)) ? aff[x][y+1][z][1] : low;
                F posz = (z<(zdim-1)) ? aff[x][y][z+1][2] : low;

                F m = std::max({negx,negy,negz,posx,posy,posz});

                // Don't follow the edge if it's less than low (T_l)
                if ( m > low )
                {
                    if ( negx == m || negx >= high ) { id |= 0x01; }
                    if ( negy == m || negy >= high ) { id |= 0x02; }
                    if ( negz == m || negz >= high ) { id |= 0x04; }
                    if ( posx == m || posx >= high ) { id |= 0x08; }
                    if ( posy == m || posy >= high ) { id |= 0x10; }
                    if ( posz == m || posz >= high ) { id |= 0x20; }
                }
            }


    const std::ptrdiff_t dir[6] = { -1, -xdim, -xdim*ydim, 1, xdim, xdim*ydim };
    // direction (-x,-y,-z, +x,+y,+z)
    const ID dirmask[6]  = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20 };
    // opposite direction (+x,+y,+z, -x,-y,-z)
    const ID idirmask[6] = { 0x08, 0x10, 0x20, 0x01, 0x02, 0x04 };

    // get plato escapers

    std::vector<std::ptrdiff_t> bfs;

    // For each vertex, for each direction..
    for ( std::ptrdiff_t idx = 0; idx < size; ++idx )
    {
        for ( std::ptrdiff_t d = 0; d < 6; ++d )
        {
            if ( seg_raw[idx] & dirmask[d] )
            {
                if ( !(seg_raw[idx+dir[d]] & idirmask[d]) )
                {
                    seg_raw[idx] |= 0x40;
                    bfs.push_back(idx);
                    d = 6; // break;
                }
            }
        }
    }

    std::size_t bfs_index = 0;

    while ( bfs_index < bfs.size() )
    {
        std::ptrdiff_t idx = bfs[bfs_index];

        ID to_set = 0;

        for ( std::ptrdiff_t d = 0; d < 6; ++d )
        {
            if ( seg_raw[idx] & dirmask[d] )
            {
                if ( seg_raw[idx+dir[d]] & idirmask[d] )
                {
                    if ( !( seg_raw[idx+dir[d]] & 0x40 ) )
                    {
                        bfs.push_back(idx+dir[d]);
                        seg_raw[idx+dir[d]] |= 0x40;
                    }
                }
                else
                {
                    to_set = dirmask[d];
                }
            }
        }
        seg_raw[idx] = to_set;
        ++bfs_index;
    }

    bfs.clear();

    //std::vector<std::size_t> counts({0});

    counts.resize(1);
    counts[0] = 0;

    ID next_id = static_cast<ID>(1);

    for ( std::ptrdiff_t idx = 0; idx < size; ++idx )
    {
        if ( seg_raw[idx] == 0 )
        {
            seg_raw[idx] |= traits::high_bit;
            ++counts[0];
        }

        if ( !( seg_raw[idx] & traits::high_bit ) && seg_raw[idx] )
        {
            bfs.push_back(idx);
            bfs_index = 0;
            seg_raw[idx] |= 0x40;

            while ( bfs_index < bfs.size() )
            {
                std::ptrdiff_t me = bfs[bfs_index];

                for ( std::ptrdiff_t d = 0; d < 6; ++d )
                {
                    if ( seg_raw[me] & dirmask[d] )
                    {
                        std::ptrdiff_t him = me + dir[d];
                        if ( seg_raw[him] & traits::high_bit )
                        {
                            counts[ seg_raw[him] & ~traits::high_bit ]
                                += bfs.size();

                            for ( auto& it: bfs )
                            {
                                seg_raw[it] = seg_raw[him];
                            }

                            bfs.clear();
                            d = 6; // break
                        }
                        else if ( !( seg_raw[him] & 0x40 ) )
                        {
                            seg_raw[him] |= 0x40;
                            bfs.push_back( him );

                        }
                    }
                }
                ++bfs_index;
            }

            if ( bfs.size() )
            {
                counts.push_back( bfs.size() );
                for ( auto& it: bfs )
                {
                    seg_raw[it] = traits::high_bit | next_id;
                }
                ++next_id;
                bfs.clear();
            }
        }
    }

    std::cout << "found: " << (next_id-1) << " components\n";

    for ( std::ptrdiff_t idx = 0; idx < size; ++idx )
    {
        seg_raw[idx] &= traits::mask;
    }

    return std::make_pair(seg_ptr, next_id-1);
}

template< typename ID, typename F >
inline region_graph_ptr<ID,F>
get_region_graph( const affinity_graph_ptr<F>& aff_ptr,
                  const volume_ptr<ID> seg_ptr,
                  std::size_t max_segid,
                  int type )
{
    /*
    Generate the region graph of edges connecting all pairs of neighboring
    segments, using a full affinity graph and map of labels for each pixel.
    The region graph is sorted in order of decreasing edge weight (see 'Output').

    Input:
      aff_ptr = affinity graph
      seg_ptr = vertex labels
      max_segid = # labels
      type = 1 just_watershed_avg, 2/3 just_watershed_bup/just_watershed
    Output:
      rg_ptr = region graph, which is a decreasingly sorted list of
        (max edge weight, label 1, label 2)
    */

    std::ptrdiff_t xdim = aff_ptr->shape()[0];
    std::ptrdiff_t ydim = aff_ptr->shape()[1];
    std::ptrdiff_t zdim = aff_ptr->shape()[2];

    volume<ID>& seg = *seg_ptr;
    affinity_graph<F> aff = *aff_ptr;

    region_graph_ptr<ID,F> rg_ptr( new region_graph<ID,F> );

    region_graph<ID,F>& rg = *rg_ptr;

    // just_watershed_bup, just_watershed
    if (type >= 2)
    {
        // edges = label-indexed vector, where elements are maps from labels to edge weights
        std::vector<std::map<ID,F>> edges(max_segid+1);

        // Iterate over all edges in the pixel affinity graph
        // Set edges[label1][label2] = max weight edge between segment1 and segment2
        for ( std::ptrdiff_t z = 0; z < zdim; ++z )
            for ( std::ptrdiff_t y = 0; y < ydim; ++y )
                for ( std::ptrdiff_t x = 0; x < xdim; ++x )
                {
                    if ( (x > 0) && seg[x][y][z] && seg[x-1][y][z] )
                    {
                        // minmax ensures `edges` contains no duplicate edges
                        auto mm = std::minmax(seg[x][y][z], seg[x-1][y][z]);
                        F& curr = edges[mm.first][mm.second];
                        curr = std::max(curr, aff[x][y][z][0]);
                    }
                    if ( (y > 0) && seg[x][y][z] && seg[x][y-1][z] )
                    {
                        auto mm = std::minmax(seg[x][y][z], seg[x][y-1][z]);
                        F& curr = edges[mm.first][mm.second];
                        curr = std::max(curr, aff[x][y][z][1]);
                    }
                    if ( (z > 0) && seg[x][y][z] && seg[x][y][z-1] )
                    {
                        auto mm = std::minmax(seg[x][y][z], seg[x][y][z-1]);
                        F& curr = edges[mm.first][mm.second];
                        curr = std::max(curr, aff[x][y][z][2]);
                    }
                }

        // Add each max inter-segment edge to the region graph
        for ( ID id1 = 1; id1 <= max_segid; ++id1 )
            for ( const auto& p: edges[id1] )
                // Add (max edge weight b/t 1-2, label1, label2)
                rg.emplace_back(p.second, id1, p.first);
    }
    // just_watershed_avg
    else
    {
        // edges = label-indexed vector, where elements are maps from labels to
        // (edge weight, # edges considered) pairs
        std::vector<std::map<ID,std::pair<F,std::size_t>>> edges(max_segid+1);

        // Iterate over all edges in the pixel affinity graph
        // Set edges[label1][label2] = (sum of edge weights, # edge weights considered)
        // for all edges between segment1 and segment2
        for ( std::ptrdiff_t z = 0; z < zdim; ++z )
            for ( std::ptrdiff_t y = 0; y < ydim; ++y )
                for ( std::ptrdiff_t x = 0; x < xdim; ++x )
                {
                    if ( (x > 0) && seg[x][y][z] && seg[x-1][y][z] )
                    {
                        // minmax ensures `edges` contains no duplicate edges
                        auto mm = std::minmax(seg[x][y][z], seg[x-1][y][z]);
                        auto& curr = edges[mm.first][mm.second];
                        curr.first += aff[x][y][z][0];
                        ++curr.second;
                    }
                    if ( (y > 0) && seg[x][y][z] && seg[x][y-1][z] )
                    {
                        auto mm = std::minmax(seg[x][y][z], seg[x][y-1][z]);
                        auto& curr = edges[mm.first][mm.second];
                        curr.first += aff[x][y][z][1];
                        ++curr.second;
                    }
                    if ( (z > 0) && seg[x][y][z] && seg[x][y][z-1] )
                    {
                        auto mm = std::minmax(seg[x][y][z], seg[x][y][z-1]);
                        auto& curr = edges[mm.first][mm.second];
                        curr.first += aff[x][y][z][0];
                        ++curr.second;
                    }
                }

        // Add each average inter-segment edge to the region graph
        for ( ID id1 = 1; id1 <= max_segid; ++id1 )
            for ( const auto& p: edges[id1] )
                // Add (avg edge weight b/t 1-2, label1, label2)
                rg.emplace_back(p.second.first/p.second.second, id1, p.first);
    }

    std::cout << "Region graph size: " << rg.size() << std::endl;

    std::sort(std::begin(rg), std::end(rg),
              std::greater<std::tuple<F,ID,ID>>());

    std::cout << "Sorted" << std::endl;

    return rg_ptr;
}


template< typename ID, typename F, typename L >
inline void yet_another_watershed( const volume_ptr<ID>& seg_ptr,
                                   const region_graph_ptr<ID,F> rg_ptr,
                                   std::vector<std::size_t>& counts,
                                   const L& lowl)
{
    /*
    Iterate through all pairs of segments in region graph, merging pairs
    with edge weight > `lowl` (T_e)
    
    Remark: Merges of groups of >2 segments are made
    when the edges connecting segments in that group have weight > `lowl`.

    Remark: This seems like a duplicate of merge_segments_with_function,
    except with `func` allowing all merges regardless of size, so long as
    f > lowl?

    Inputs:
      seg_ptr = labels of each vertex
      rg_ptr = region graph, containing (edge weight, label 1, label 2)
      counts = size of each segment
      lowl = T_e
    Outputs:
      None, but seg_ptr and rg_ptr are updated.
    */

    F low = static_cast<F>(lowl);

    // new_counts = new sizes of segments
    // remaps = new labels of segments
    std::vector<std::size_t> new_counts({0});
    std::vector<ID>          remaps(counts.size());

    // sets = Quick-union structure
    // maxs = New edges?
    zi::disjoint_sets<ID>    sets(counts.size());
    std::vector<F>           maxs(counts.size());


    region_graph<ID,F>& rg  = *rg_ptr;

    // next_id = Next label
    // merged = Count of # segments merged
    ID next_id = 1;
    ID merged = 0;

    //
    // Merge all pairs of segments in `rg` with edge weight > `low` (T_e).

    for ( auto& it: rg )
    {
        // The T_e check.
        if ( std::get<0>(it) <= low )
            break;

        // s1, s2 = labels of the two segments connected by this edge
        ID s1 = std::get<1>(it);
        ID s2 = std::get<2>(it);

        // f = weight of this edge
        F f = std::get<0>(it);

        // Summary: If the two segments are both part of a valid region,
        // then merge the two segments.
        if ( s1 && s2 )
        {
            if ( (remaps[s1] == 0) || (remaps[s2] == 0) )
            {
                if ( remaps[s1] == 0 )
                {
                    std::swap(s1,s2);
                }

                // If both segments haven't been part of a merge yet
                if ( remaps[s1] == 0 )
                {
                    maxs[next_id] = f;
                    remaps[s1] = remaps[s2] = next_id;
                    new_counts.push_back(counts[s1]+counts[s2]);
                    ++next_id;
                }
                // If only 1 segment has been part of a merge so far
                else
                {
                    ID actual = sets.find_set(remaps[s1]);
                    remaps[s2] = remaps[s1];
                    new_counts[actual] += counts[s2];
                }
            }
            // If both segments have been part of a merge already
            else
            {
                ID a1 = sets.find_set(remaps[s1]);
                ID a2 = sets.find_set(remaps[s2]);

                if ( 0 && a1 != a2 && ((maxs[a1]==f)||(maxs[a2]==f)) )
                {
                    ++merged;
                    new_counts[a1] += new_counts[a2];
                    new_counts[a2] = 0;
                    maxs[a1] = std::max(maxs[a1],maxs[a2]);
                    maxs[a2] = 0;
                    ID a = sets.join(a1,a2);
                    std::swap(new_counts[a], new_counts[a1]);
                    std::swap(maxs[a], maxs[a1]);
                }
            }
        }
    }

    //
    // Compress the representation of the above merging (i.e., get rid of the
    // quick union stuff)

    std::vector<ID> remaps2(counts.size());

    next_id = 1;

    for ( ID id = 0; id < counts.size(); ++id )
    {
        ID s = sets.find_set(remaps[id]);
        if ( s && (remaps2[s]==0) )
        {
            remaps2[s] = next_id;
            new_counts[next_id] = new_counts[s];
            ++next_id;
        }
    }

    new_counts.resize(next_id);


    //
    // Update seg_ptr with new (merged) labels of segments

    std::ptrdiff_t xdim = seg_ptr->shape()[0];
    std::ptrdiff_t ydim = seg_ptr->shape()[1];
    std::ptrdiff_t zdim = seg_ptr->shape()[2];

    std::ptrdiff_t total = xdim * ydim * zdim;

    ID* seg_raw = seg_ptr->data();

    for ( std::ptrdiff_t idx = 0; idx < total; ++idx )
    {
        seg_raw[idx] = remaps2[remaps[seg_raw[idx]]];
    }

    //
    // Update rg_ptr with new edge weights

    xdim = seg_ptr->shape()[0];
    region_graph<ID,F> new_rg;

    for ( auto& it: rg )
    {
        ID s1 = remaps2[remaps[std::get<1>(it)]];
        ID s2 = remaps2[remaps[std::get<2>(it)]];

        if ( s1 != s2 && s1 && s2 )
        {
            auto mm = std::minmax(s1,s2);
            new_rg.emplace_back(std::get<0>(it), mm.first, mm.second);
        }
    }

    rg.swap(new_rg);

    counts.swap(new_counts);

    std::cout << "Done with region watershed, total: " << counts.size() << std::endl;

    std::cout << "Done with updating the region graph, size: "
              << rg.size() << std::endl;
}


template< typename ID, typename F, typename FN, typename M >
inline void merge_segments_with_function( const volume_ptr<ID>& seg_ptr,
                                          const region_graph_ptr<ID,F> rg_ptr,
                                          std::vector<std::size_t>& counts,
                                          const FN& func,
//                                          const L& tholds,
                                          const M& lowt )
{
    /*
    Iterate through all pairs of segments in region graph, merging pairs
    with at least one segment's size < func( edge weight b/t segments ).
    
    Remark: Merges of groups of >2 segments are made
    according to sizes of segments that are updated as merges are being made.
    (Otherwise, sizes would be underreported.)

    Inputs:
      seg_ptr = labels of each vertex
      rg_ptr = region graph, containing (edge weight, label 1, label 2)
      counts = size of each segment
      lowt = min size for segments, after merge
    Outputs:
      None, but seg_ptr and rg_ptr are updated.
    */

    // Quick-union structure
    zi::disjoint_sets<ID> sets(counts.size());

    region_graph<ID,F>& rg  = *rg_ptr;

    // Iterate through pairs in the region graph..
    // Merge all pairs with at least one size 
    for ( auto& it: rg )
    {
        // max size for at least one segment = func( edge weight b/t segments )
        std::size_t size = func(std::get<0>(it));

        if ( size == 0 )
            break;

        ID s1 = sets.find_set(std::get<1>(it));
        ID s2 = sets.find_set(std::get<2>(it));

        // Merge if either segment's size < 'size'
        if ( s1 != s2 && s1 && s2 )
            if ( (counts[s1] < size) || (counts[s2] < size) )
            {
                counts[s1] += counts[s2];
                counts[s2]  = 0;
                ID s = sets.join(s1,s2);
                std::swap(counts[s], counts[s1]);
            }
    }

    // This was from merge_segments(..., tholds, ...)
    // The uncommented version is merge_segments_with_function(..., func, ...)
    /*
    typename region_graph<ID,F>::iterator rit = rg_ptr->begin();

    for ( auto& it: tholds )
    {
        std::size_t size = static_cast<std::size_t>(it.first);
        F           thld = static_cast<F>(it.second);

        while ( (rit != rg.end()) && ( std::get<0>(*rit) > thld) )
        {
            ID s1 = sets.find_set(std::get<1>(*rit));
            ID s2 = sets.find_set(std::get<2>(*rit));

            if ( s1 != s2 && s1 && s2 )
            {
                if ( (counts[s1] < size) || (counts[s2] < size) )
                {
                    counts[s1] += counts[s2];
                    counts[s2]  = 0;
                    ID s = sets.join(s1,s2);
                    std::swap(counts[s], counts[s1]);
                }
            }
            ++rit;
        }
    }
    */

    std::cout << "Done with merging" << std::endl;


    //
    // Update the segment graph of pixel labels

    std::vector<ID> remaps(counts.size());

    ID next_id = 1;

    std::size_t low = static_cast<std::size_t>(lowt);

    // Discard all elements with size smaller than `low`
    for ( ID id = 0; id < counts.size(); ++id )
    {
        ID s = sets.find_set(id);
        if ( s && (remaps[s] == 0) && (counts[s] >= low) )
        {
            remaps[s] = next_id;
            counts[next_id] = counts[s];
            ++next_id;
        }
    }

    counts.resize(next_id);

    std::ptrdiff_t xdim = seg_ptr->shape()[0];
    std::ptrdiff_t ydim = seg_ptr->shape()[1];
    std::ptrdiff_t zdim = seg_ptr->shape()[2];

    std::ptrdiff_t total = xdim * ydim * zdim;

    ID* seg_raw = seg_ptr->data();

    for ( std::ptrdiff_t idx = 0; idx < total; ++idx )
    {
        seg_raw[idx] = remaps[sets.find_set(seg_raw[idx])];
    }

    std::cout << "Done with remapping, total: " << (next_id-1) << std::endl;


    //
    // Update the region graph of edges

    region_graph<ID,F> new_rg;

    for ( auto& it: rg )
    {
        ID s1 = remaps[sets.find_set(std::get<1>(it))];
        ID s2 = remaps[sets.find_set(std::get<2>(it))];

        if ( s1 != s2 && s1 && s2 )
        {
            auto mm = std::minmax(s1,s2);
            new_rg.emplace_back(std::get<0>(it), mm.first, mm.second);
        }
    }

    rg.swap(new_rg);

    std::cout << "Done with updating the region graph, size: "
              << rg.size() << std::endl;
}


template< typename T >
affinity_graph_ptr<T> mult_aff( const affinity_graph_ptr<T>& aff, int n )
{

    std::ptrdiff_t xdim = aff->shape()[0];
    std::ptrdiff_t ydim = aff->shape()[1];
    std::ptrdiff_t zdim = aff->shape()[2];

    affinity_graph_ptr<T> r(new affinity_graph<T>
                            (boost::extents[xdim*n][ydim*n][zdim*n][3],
                             boost::fortran_storage_order()));

    for ( std::ptrdiff_t z = 0; z < n; ++z )
        for ( std::ptrdiff_t y = 0; y < n; ++y )
            for ( std::ptrdiff_t x = 0; x < n; ++x )
            {

                typedef boost::multi_array_types::index_range range;

                (*r)[boost::indices[range(x*xdim,(x+1)*xdim)]
                     [range(y*ydim,(y+1)*ydim)]
                     [range(z*zdim,(z+1)*zdim)]
                     [range(0,3)]] = *aff;

            }

    return r;

}


/*
Input:
  v = edge weight between pair of segments
Output:
  max size (voxels) of at least one segment, in order for that pair of
  segments to be merged.
*/

std::size_t limit_fn_avg( float v )
{
    // size threshold based on affinity
    if ( v > 1 )
        printf("Bad v: %f", v);
        return 2000;

    if ( v < 0.5 )
        return 25;

    /* XXX: Uh.. this isn't smooth? */
    v *= 7;

    return static_cast<std::size_t>(v*v*v*v);
}

std::size_t limit_fn_bup( float v )
{
    // size threshold based on affinity
    if ( v > 1 )
        printf("Bad v: %f", v);
        return 2000;

    if ( v < 0.5 )
        return 25;

    /* XXX: Uh.. this isn't smooth? */
    v *= 4;

    return static_cast<std::size_t>(v*v*v*v*v);
}


int main( int argc, char *argv[] )
{
    /*
    Accepts argument `infile`.
        Path to input affinity graph.
    Accepts argument `infile_size`.
        Dimensions of infile, as 1 side of a cube.
    Accepts argument `type`.
        ALL: CC w/ T_h=0.99, then WS w/ T_l=0.3
        0: (nothing else)
        1: just_watershed_avg (+merge pairs w/ size > limit_fn_avg(avg edge))
        2: just_watershed_bup (+merge pairs w/ size > limig_fn_bup(max edge))
        3: just_watershed (+merge pairs w/ max edge > 0.1)
    */

    if (argc < 4)
        throw std::invalid_argument( "insufficient arguments" );

    // Detect input file.
    char* infile = argv[1];
    int infile_size = std::stoi(argv[2]);
    // Detect type.
    int type = std::stoi(argv[3]);
    if (type > 3 || type < 0)
        throw std::invalid_argument( "invalid argument" );


    auto start = std::chrono::system_clock::now();
    auto start_master = start;
    typedef std::chrono::duration<int> seconds_type;


    //std::cout << std::hex << watershed_traits<uint32_t>::high_bit << std::endl;

    affinity_graph_ptr<float> aff;
    aff = read_affinity_graph_from_file<float>(infile, infile_size, infile_size, infile_size);

    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start).count() << "ms: "
            << "data loaded" << std::endl;
    start = std::chrono::system_clock::now();


    //std::fill_n(aff->data(), 160*160*160*3, 0);
    // (*aff)[10][10][10][0] = 0.5;
    // (*aff)[10][10][150][0] = 0.5;
    // (*aff)[10][100][100][0] = 0.5;
    //aff = mult_aff(aff, 3);
    //std::cout << "Multiplied" << std::endl;


    std::vector<std::size_t> counts;

    // seg is a pair, containing (vertex labels, # labels)
    // Cool! same return value as my watershed function lol
    auto seg = simple_watershed<uint32_t>(aff, 0.3, 0.99, counts);


    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start).count() << "ms: "
            << "simple ws done" << std::endl;
    start = std::chrono::system_clock::now();


    if (type > 0)
    {

    // Pass (affinity graph, vertex labels, # labels) into get_region_graph
    // Receive a region graph that is a list of (weight of max edge, label 1, label 2)
    // in decreasing order of edge weight.
    auto rg = get_region_graph<uint32_t,float>(aff, seg.first, seg.second, type);


    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start).count() << "ms: "
            << "region graph created" << std::endl;
    start = std::chrono::system_clock::now();


    if (type == 1)
    {
        merge_segments_with_function(seg.first, rg, counts, limit_fn_avg, 25);
    }
    else if (type == 2)
    {
        merge_segments_with_function(seg.first, rg, counts, limit_fn_bup, 25);
    }
    else if (type == 3)
    {
        // Note that the segment size threshold T_s isn't implemented here.
        yet_another_watershed(seg.first, rg, counts, 0.1);
    }

    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start).count() << "ms: "
            << "region segments watersheded/merged" << std::endl;

    }


    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start_master).count() << "ms total: "
            << "finished, now outputting..." << std::endl;

    char outfile[100];
    sprintf(outfile, "aleks-%d-%d.out", infile_size, type);
    write_volume_to_file(outfile, seg.first);

}
