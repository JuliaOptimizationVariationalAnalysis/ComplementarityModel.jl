for field in [:Cmin, :CFB]
  fin = Symbol(field, "!")
  @eval begin
    @doc """
        $($field)(Fx, Gx)
    Compute the C-function value.
    """
    function $field(Fx::T, Gx::T) where T
      Cx = similar(Fx)
      return $fin(Cx, Fx, Gx)
    end
  end
  @eval export $field, $fin
end

function Cmin!(Cx, Fx, Gx)
  Cx .= min.(Fx, Gx)
  return Cx
end

function CFB!(Cx, Fx, Gx) # Fis
  @. Cx = Fx + Gx - sqrt(Fx^2 + Gx^2) 
  return Cx
end