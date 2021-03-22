package org.pankratzlab.supernovo.output;

import java.lang.reflect.Field;
import java.util.stream.Collector;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import org.pankratzlab.supernovo.App;
import com.google.common.base.Optional;

/**
 * Interface to specify and allow a class's public fields to be used as output columns instead of
 * resolving to a single column with {@link Object#toString()}
 */
public interface OutputFields {

  /** Holds static final fields that should not be printed on output */
  static class Constants {
    private Constants() {}

    static final String DELIM = "\t";
    static final String MISSING = ".";
    static final Collector<CharSequence, ?, String> JOIN_COLLECTOR = Collectors.joining(DELIM);
  }

  default String generateLine() {
    return fieldValues().collect(Constants.JOIN_COLLECTOR);
  }

  default Stream<String> fieldValues() {
    return Stream.of(this.getClass().getFields())
        .map(this::getOwnField)
        .flatMap(this::recurseValues);
  }

  default Object getOwnField(Field field) {
    try {
      return field.get(this);
    } catch (IllegalArgumentException | IllegalAccessException e) {
      throw new IllegalStateException(e);
    }
  }

  default Stream<String> recurseValues(Object value) {
    if (value instanceof OutputFields) {
      return ((OutputFields) value).fieldValues();
    }
    if (value instanceof Optional<?>) {
      return Stream.of(((Optional<?>) value).transform(Object::toString).or(Constants.MISSING));
    }
    if (value == null) {
      return Stream.of(Constants.MISSING);
    }
    return Stream.of(value.toString());
  }

  default String generateHeader() {
    return fieldHeaders().collect(Constants.JOIN_COLLECTOR);
  }

  default Stream<String> fieldHeaders() {
    return Stream.of(this.getClass().getFields()).flatMap(this::recurseHeaders);
  }

  default Stream<String> recurseHeaders(Field field) {
    Class<?> fieldType = field.getType();
    if (OutputFields.class.isAssignableFrom(fieldType)) {
      try {
        return prefixHeaders(field, ((OutputFields) field.get(this)).fieldHeaders());
      } catch (IllegalArgumentException | IllegalAccessException e) {
        App.LOG.error(e);
      }
    }

    return Stream.of(field.getName());
  }

  static Stream<String> prefixHeaders(Field field, Stream<? extends Object> headers) {
    final String prefix = field.getName() + "_";
    return headers.map(Object::toString).map(h -> prefix + h);
  }
}
