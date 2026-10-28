use serde::{Deserialize, Deserializer, Serialize, Serializer};

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Number {
    Finite(f64),
    Na,
    Nan,
    PosInf,
    NegInf,
}

impl Number {
    pub fn from_value(value: f64) -> Self {
        if value.is_nan() {
            Self::Nan
        } else if value == f64::INFINITY {
            Self::PosInf
        } else if value == f64::NEG_INFINITY {
            Self::NegInf
        } else {
            Self::Finite(value)
        }
    }

    pub fn value(self) -> f64 {
        match self {
            Self::Finite(value) => value,
            Self::Na | Self::Nan => f64::NAN,
            Self::PosInf => f64::INFINITY,
            Self::NegInf => f64::NEG_INFINITY,
        }
    }

    pub fn text(self) -> String {
        match self {
            Self::Finite(value) => value.to_string(),
            Self::Na => "NA".into(),
            Self::Nan => "NaN".into(),
            Self::PosInf => "Inf".into(),
            Self::NegInf => "-Inf".into(),
        }
    }
}

impl Serialize for Number {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        match self {
            Self::Finite(value) if value.is_finite() => serializer.serialize_f64(*value),
            Self::Finite(_) => Err(serde::ser::Error::custom("nonfinite Finite state")),
            _ => serializer.serialize_str(&self.text()),
        }
    }
}

impl<'de> Deserialize<'de> for Number {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        #[derive(Deserialize)]
        #[serde(untagged)]
        enum Encoded {
            Value(f64),
            Token(String),
        }
        match Encoded::deserialize(deserializer)? {
            Encoded::Value(value) if value.is_finite() => Ok(Self::Finite(value)),
            Encoded::Token(token) => match token.as_str() {
                "NA" => Ok(Self::Na),
                "NaN" => Ok(Self::Nan),
                "Inf" => Ok(Self::PosInf),
                "-Inf" => Ok(Self::NegInf),
                _ => Err(serde::de::Error::custom("invalid numeric-state token")),
            },
            _ => Err(serde::de::Error::custom("invalid finite numeric state")),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "shape", content = "values", rename_all = "kebab-case")]
pub enum RowValues {
    Scalar(f64),
    Rows(Vec<f64>),
}
