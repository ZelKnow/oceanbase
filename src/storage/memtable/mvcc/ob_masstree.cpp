
#include "ob_masstree.h"
#include "lib/allocator/ob_qsync.h"
#include "common/object/ob_object.h"
#include "lib/allocator/ob_retire_station.h"
#include "lib/oblog/ob_log_module.h"
#include "share/ob_errno.h"


namespace oceanbase
{
namespace keybtree
{
using namespace oceanbase::common;
using namespace Masstree;

namespace masstree
{
template<typename MassTreeKey, typename MassTreeVal>
int ObMassTree<MassTreeKey, MassTreeVal>::insert(const MassTreeKey key, MassTreeVal &value) {
    const ObObj *obj = key.get_ptr();
    if (obj->get_type() == ObObjType::ObIntType) {
        cursor_type 
    } else if (obj->get_type() == ObObjType::ObCharType) {

    }
}
}
}
}