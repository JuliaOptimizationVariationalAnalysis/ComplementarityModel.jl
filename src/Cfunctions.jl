for field in [:Cmin, :CFB]
  fin = Symbol(field, "!")
  @eval begin
    @doc """
        $($field)(Fx, Gx)
    Compute the C-function value.
    """
    function $field(Fx, Gx)
      Cx = similar(Fx)
      return $fin(Cx, Fx, Gx)
    end
  end
  Jf = Symbol("J", field)
  Jfin = Symbol("J", field, "!")
  @eval begin
    @doc """
        $($Jf)(Fx, Gx, JFx, JGx)
    Compute the jacobian's diagonal of the C-function.
    """
    function $Jf(Fx, Gx, JFx, JGx)
      Cx = similar(JFx)
      return $Jfin(Cx, Fx, Gx, JFx, JGx)
    end
  end
  Jfinv = Symbol("J", field, "v!")
  Jfv = Symbol("J", field, "v")
  @eval begin
    @doc """
        $($Jfv)(Fx, Gx, JFx, JGx)
    Compute the jacobian's diagonal of the C-function.
    """
    function $Jfv(Fx, Gx, JFx, JGx, v)
      Cx = similar(Fx)
      return $Jfinv(Cx, Fx, Gx, JFx, JGx, v)
    end
  end
  @eval export $field, $fin, $Jf, $Jfin, $Jfv, $Jfinv
end

function Cmin!(Cx, Fx, Gx)
  Cx .= min.(Fx, Gx)
  return Cx
end

function JCmin!(Cx, Fx, Gx, JFx, JGx)
  for i=1:size(Cx, 1)
    @views Cx[i, :] .= Fx[i] ≤ Gx[i] ? JFx[i, :] : JGx[i, :]
  end
  return Cx
end

function JCminv!(w, Fx, Gx, JFx, JGx, v)
  for i=1:length(w)
    w[i] = Fx[i] ≤ Gx[i] ? dot(view(JFx, i, :), v) : dot(view(JGx, i, :), v)
  end
  return w
end

function CFB!(Cx, Fx, Gx) # Fischer-Burmeister C-function
  @. Cx = Fx + Gx - sqrt(Fx^2 + Gx^2) 
  return Cx
end

function JCFB!(Cx, Fx, Gx, JFx, JGx)
  for i=1:size(Cx, 1)
    if (Fx[i].^2 .+ Gx[i].^2) == 0
      view(Cx, i, :) .= view(JFx, i, :) .+ view(JGx, i, :) .- (view(JFx, i, :) .* Fx[i] + view(JGx, i, :) .* Gx[i]) ./ sqrt(Fx[i].^2 .+ Gx[i].^2)
    else
      view(Cx, i, :) .= view(JFx, i, :) .+ view(JGx, i, :) # the jacobian when F=G=0 is not correct
    end
  end
  return Cx
end

function JCFBv!(w, Fx, Gx, JFx, JGx, v)
  for i=1:length(w)
    if (Fx[i].^2 .+ Gx[i].^2) == 0
        w[i] = dot(view(JFx, i, :), v) + dot(view(JGx, i, :), v) - (dot(view(JFx, i, :), v) * Fx[i] + dot(view(JGx, i, :), v) * Gx[i]) / sqrt(Fx[i]^2 .+ Gx[i]^2)
    else
        w[i] = dot(view(JFx, i, :), v) + dot(view(JGx, i, :), v) # the jacobian when F=G=0 is not correct
    end
  end
  return w
end
