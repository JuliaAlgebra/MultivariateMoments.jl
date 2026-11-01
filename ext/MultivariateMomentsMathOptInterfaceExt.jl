module MultivariateMomentsMathOptInterfaceExt

import MathOptInterface as MOI
import MultivariateMoments as MM

const ObjectWithoutIndex = Union{
    MM.AbstractMomentMatrix{<:MOI.Utilities.ObjectWithoutIndex},
    MM.AbstractMeasure{<:MOI.Utilities.ObjectWithoutIndex},
}
# Empty tuples are already handled by MOI, independently of this extension.
const ObjectOrTupleWithoutIndex = Union{
    ObjectWithoutIndex,
    Tuple{ObjectWithoutIndex,Vararg{ObjectWithoutIndex}},
}
const ObjectOrTupleOrArrayWithoutIndex =
    Union{ObjectOrTupleWithoutIndex,AbstractArray{<:ObjectOrTupleWithoutIndex}}

function MOI.Utilities.map_indices(
    ::F,
    x::ObjectOrTupleOrArrayWithoutIndex,
) where {F<:Function}
    return x
end

function MOI.Utilities.substitute_variables(
    ::F,
    x::ObjectOrTupleOrArrayWithoutIndex,
) where {F<:Function}
    return x
end

end
